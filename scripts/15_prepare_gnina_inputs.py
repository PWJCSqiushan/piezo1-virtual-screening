from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.docking import (
    docking_box_from_atoms,
    read_candidate_list,
    read_ranked_compounds,
    receptor_policy_summary,
    write_receptor_pdb,
)
from piezo_vs.drugclip_io import extract_pocket_atoms, parse_residue_tokens, select_consensus_row
from piezo_vs.io_utils import read_json, sha256_file, write_json


CANDIDATE_MANIFEST_FIELDS = [
    "candidate_index",
    "compound_id",
    "canonical_smiles",
    "site_direction",
    "pdb_id",
    "consensus_id",
    "source",
    "source_id",
    "source_url",
    "rank",
    "drugclip_score",
    "review_status",
    "approved",
    "evaluation_complete",
    "preparation_status",
    "failure_reason", "assessment_evidence_sha256", "reviewer", "reviewed_at_utc", "assessment_comment", "selection_reason", "source_ids", "source_tasks",
]


def prepare_ligand(smiles: str, seed: int):
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:
        raise RuntimeError("GNINA input preparation requires RDKit") from exc
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("RDKit rejected the ranked SMILES")
    molecule = Chem.AddHs(molecule)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = seed
    parameters.numThreads = 1
    if AllChem.EmbedMolecule(molecule, parameters) != 0:
        raise ValueError("RDKit could not generate a 3D conformer")
    if AllChem.MMFFHasAllMoleculeParams(molecule):
        AllChem.MMFFOptimizeMolecule(molecule, maxIters=500)
    else:
        AllChem.UFFOptimizeMolecule(molecule, maxIters=500)
    return molecule


def _candidate_rows(args, pdb_id: str, consensus_id: str) -> tuple[list[dict[str, str]], dict[str, object]]:
    if not args.ranked_csv and not args.candidate_list:
        raise ValueError("Provide --candidate-list for explicit candidates or --ranked-csv for the legacy demo")
    ranked_rows: list[dict[str, str]] | None = None
    if args.ranked_csv:
        ranked_rows = read_ranked_compounds(
            args.ranked_csv,
            pdb_id,
            consensus_id,
            top_n=None if args.candidate_list else args.top_n,
        )
    if args.candidate_list:
        if args.purpose == "formal_screening":
            selection_path = args.candidate_list.parent / "selection_manifest.json"
            selection = read_json(selection_path)
            if selection.get("selected_files",{}).get(args.candidate_list.name) != sha256_file(args.candidate_list):
                raise ValueError("Candidate list does not match the selection manifest")
            if selection["assessment_sha256"] != sha256_file(args.candidate_list.parent / "candidate_assessment.csv"):
                raise ValueError("Selection assessment snapshot changed")
        rows = read_candidate_list(
            args.candidate_list,
            pdb_id,
            consensus_id,
            purpose=args.purpose,
            ranked_rows=ranked_rows,
        )
        source = {
            "type": "explicit_candidate_list",
            "path": str(args.candidate_list.resolve()),
            "sha256": sha256_file(args.candidate_list),
        }
        if args.ranked_csv:
            source["ranked_csv_path"] = str(args.ranked_csv.resolve())
            source["ranked_csv_sha256"] = sha256_file(args.ranked_csv)
        return rows, source
    if args.purpose == "formal_screening":
        raise ValueError("Formal GNINA requires --candidate-list; --ranked-csv --top-n is legacy technical mode")
    assert ranked_rows is not None
    rows = []
    for row in ranked_rows:
        row = dict(row)
        row["canonical_smiles"] = row.get("canonical_smiles") or row.get("smiles", "")
        row.setdefault("site_direction", args.site_direction or consensus_id)
        row.setdefault("review_status", "legacy_ranked")
        row.setdefault("approved", "")
        row.setdefault("evaluation_complete", "")
        rows.append(row)
    return rows, {
        "type": "legacy_ranked_csv",
        "path": str(args.ranked_csv.resolve()),
        "sha256": sha256_file(args.ranked_csv),
        "requested_top_n": args.top_n,
    }


def _write_rows(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare traceable GNINA receptor, explicit candidate ligands and search box."
    )
    parser.add_argument("--pdb-id", required=True)
    parser.add_argument("--consensus-id", required=True)
    parser.add_argument("--ranked-csv", type=Path)
    parser.add_argument(
        "--candidate-list",
        type=Path,
        help=(
            "Explicit CSV with compound_id, canonical_smiles, site_direction, pdb_id, "
            "consensus_id, review_status and evaluation_complete; no implicit top-N is applied."
        ),
    )
    parser.add_argument("--purpose", choices=("technical_validation", "formal_screening"), default="technical_validation")
    parser.add_argument("--top-n", type=int, default=20, help="Legacy ranked CSV limit; ignored for --candidate-list")
    parser.add_argument("--site-direction", default="", help="Site label for legacy ranked CSV rows")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--allow-non-piezo1-scope", action="store_true")
    parser.add_argument("--allow-large-box-clipping", action="store_true")
    parser.add_argument(
        "--heteroatom-policy",
        choices=("remove_all", "keep_water_ions", "keep_all"),
        default="remove_all",
        help="HETATM handling is independent of pocket extraction and is recorded in the manifest.",
    )
    parser.add_argument(
        "--receptor-review-approved",
        "--approve-receptor-environment",
        dest="receptor_review_approved",
        action="store_true",
        help="Record manual approval of protein/HETATM environment for a formal run.",
    )
    parser.add_argument("--receptor-review-note", default="")
    parser.add_argument("--residue-policy", choices=("two_plus", "all_supporting_tools"), default="two_plus")
    args = parser.parse_args()
    if not args.candidate_list and not 1 <= args.top_n <= 100:
        raise SystemExit("top-n must be between 1 and 100 for legacy ranked CSV mode")
    if args.candidate_list and args.top_n < 1:
        raise SystemExit("top-n must be positive even when --candidate-list is supplied")

    pdb_id = args.pdb_id.upper()
    consensus_id = args.consensus_id.upper()
    output_dir = args.output_dir.resolve()

    structure_config = read_json(ROOT / "config" / "structures.json")
    try:
        structure = next(item for item in structure_config["structures"] if item["pdb_id"] == pdb_id)
    except StopIteration as exc:
        raise ValueError(f"{pdb_id} is not an active configured structure") from exc
    consensus_path = ROOT / "results" / "consensus" / f"{pdb_id}_consensus_pockets.csv"
    source_pdb = ROOT / "data" / "raw" / "structures" / f"{pdb_id}.pdb"
    for path in (consensus_path, source_pdb):
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.ranked_csv and not args.ranked_csv.is_file():
        raise FileNotFoundError(args.ranked_csv)
    if args.candidate_list and not args.candidate_list.is_file():
        raise FileNotFoundError(args.candidate_list)

    consensus_rules = read_json(ROOT / "config" / "consensus_rules.json")
    tier_by_support = {
        int(definition["support_count"]): tier
        for tier, definition in consensus_rules["tiers"].items()
    }
    consensus = select_consensus_row(
        consensus_path,
        pdb_id,
        consensus_id,
        tier_by_support=tier_by_support,
        required_tools=set(consensus_rules["required_tools"]),
    )
    scope = consensus.get("pocket_scope", "")
    if scope != "piezo1" and not args.allow_non_piezo1_scope:
        raise ValueError(
            f"Pocket scope {scope!r} requires explicit expert approval and --allow-non-piezo1-scope"
        )
    if args.purpose == "formal_screening" and not args.receptor_review_approved:
        raise ValueError(
            "Formal GNINA requires --receptor-review-approved after the protein/HETATM environment "
            "has been reviewed"
        )

    if args.purpose == "formal_screening" and not args.receptor_review_note.strip():
        raise ValueError("Formal receptor review needs a reviewer/date and rationale in --receptor-review-note")
    residue_field = {"two_plus": "core_residues_2plus", "all_supporting_tools": "core_residues_all_supporting_tools"}[args.residue_policy]
    residues = parse_residue_tokens(consensus[residue_field])
    pocket_atoms, missing = extract_pocket_atoms(source_pdb, residues)
    if missing:
        raise ValueError(f"{len(missing)} consensus residues are absent from the source PDB")
    center = tuple(float(consensus[key]) for key in ("center_x", "center_y", "center_z"))
    box = docking_box_from_atoms(pocket_atoms, center)
    if box["clipped_to_maximum"] and not args.allow_large_box_clipping:
        raise ValueError(
            "The pocket-derived box exceeds 40 angstrom on at least one axis. "
            "Review the pocket and use --allow-large-box-clipping only if the capped box is intended."
        )

    rows, candidate_source = _candidate_rows(args, pdb_id, consensus_id)
    if not rows:
        raise ValueError("No GNINA candidates were selected")
    for row in rows:
        row["pdb_id"] = pdb_id
        row["consensus_id"] = consensus_id
        row["canonical_smiles"] = row.get("canonical_smiles") or row.get("smiles", "")

    # Keep the complete configured protein environment.  Pocket atoms only
    # control the search box and never remove a protein chain from the receptor.
    receptor_chains = set(structure["auth_piezo1_chains"]) | set(structure.get("auth_mdfic_chains", []))
    receptor_summary = receptor_policy_summary(
        source_pdb,
        receptor_chains,
        heteroatom_policy=args.heteroatom_policy,
    )
    receptor_policy = {
        **receptor_summary,
        "environment": "all_configured_piezo1_and_mdfic_protein_chains",
        "chains": sorted(receptor_chains),
        "review_required_for_formal": True,
        "review_status": "approved" if args.receptor_review_approved else "pending",
        "review_note": args.receptor_review_note,
    }

    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    receptor_path = output_dir / "receptor.pdb"
    receptor_atom_count = write_receptor_pdb(
        source_pdb,
        receptor_path,
        receptor_chains,
        heteroatom_policy=args.heteroatom_policy,
    )

    try:
        from rdkit import Chem, rdBase
    except ImportError as exc:
        raise RuntimeError("GNINA input preparation requires RDKit") from exc
    ligand_path = output_dir / "ligands.sdf"
    writer = Chem.SDWriter(str(ligand_path))
    accepted: list[dict[str, str]] = []
    rejected: list[dict[str, str]] = []
    manifest_rows: list[dict[str, str]] = []
    try:
        for index, original_row in enumerate(rows, 1):
            row = dict(original_row)
            row["candidate_index"] = str(index)
            row["preparation_status"] = "failed"
            row["failure_reason"] = ""
            try:
                if not row.get("canonical_smiles"):
                    raise ValueError("candidate has no canonical_smiles")
                molecule = prepare_ligand(row["canonical_smiles"], args.seed + index - 1)
                molecule.SetProp("_Name", row["compound_id"])
                for key in (
                    "candidate_index",
                    "compound_id",
                    "canonical_smiles",
                    "site_direction",
                    "pdb_id",
                    "consensus_id",
                    "source",
                    "source_id",
                    "source_url",
                    "rank",
                    "drugclip_score",
                    "review_status",
                    "approved",
                    "evaluation_complete",
                ):
                    molecule.SetProp(key, row.get(key, ""))
                writer.write(molecule)
                row["preparation_status"] = "prepared"
                accepted.append(row)
            except (ValueError, RuntimeError, OverflowError) as exc:
                row["failure_reason"] = str(exc)
                rejected.append(
                    {
                        "compound_id": row.get("compound_id", ""),
                        "canonical_smiles": row.get("canonical_smiles", ""),
                        "reason": str(exc),
                    }
                )
            manifest_rows.append({field: row.get(field, "") for field in CANDIDATE_MANIFEST_FIELDS})
    finally:
        writer.close()

    rejected_path = output_dir / "rejected_ligands.csv"
    _write_rows(rejected_path, rejected, ["compound_id", "canonical_smiles", "reason"])
    candidate_manifest_path = output_dir / "candidate_manifest.csv"
    _write_rows(candidate_manifest_path, manifest_rows, CANDIDATE_MANIFEST_FIELDS)

    requested_count = len(rows)
    formal_gate_passed = all(
        row.get("approved") == "true"
        and row.get("evaluation_complete") == "true"
        and row.get("review_status") in {"approved", "accepted", "manual_approved"}
        for row in rows
    )
    technical_gate_required = args.purpose == "formal_screening" or requested_count >= 20
    config = {
        "schema_version": 2,
        "run_id": output_dir.name,
        "purpose": args.purpose,
        "pdb_id": pdb_id,
        "consensus_id": consensus_id,
        "tier": consensus["tier"],
        "support_count": int(consensus["support_count"]),
        "pocket_scope": scope,
        "residue_source_field": residue_field,
        "residue_policy": args.residue_policy,
        "box": box,
        "receptor": "receptor.pdb",
        "receptor_environment": "full_configured_protein",
        "receptor_chains": sorted(receptor_chains),
        "receptor_atom_count": receptor_atom_count,
        "receptor_policy": receptor_policy,
        "candidate_source": candidate_source,
        "candidate_manifest": "candidate_manifest.csv",
        "requested_top_n": None if args.candidate_list else args.top_n,
        "requested_candidate_count": requested_count,
        "prepared_ligand_count": len(accepted),
        "rejected_ligand_count": len(rejected),
        "preparation_status": "success" if not rejected else "partial" if accepted else "failed",
        "candidate_ids": [row["compound_id"] for row in rows],
        "seed": args.seed,
        "rdkit_version": rdBase.rdkitVersion,
        "formal_candidate_gate": {
            "required": args.purpose == "formal_screening",
            "all_evaluation_complete": all(row.get("evaluation_complete") == "true" for row in rows),
            "all_approved": all(
                row.get("approved") == "true"
                and row.get("review_status") in {"approved", "accepted", "manual_approved"}
                for row in rows
            ),
            "status": "passed" if formal_gate_passed else "not_applicable" if args.purpose != "formal_screening" else "failed",
        },
        "technical_check_gate": {
            "required_before_batch": technical_gate_required,
            "minimum_compounds": 3,
            "status": "pending" if technical_gate_required else "not_required",
            "evidence": "",
        },
        "scientific_status": "computational_docking_input_not_binding_evidence",
    }
    write_json(output_dir / "docking_config.json", config)
    output_hashes = {
        "receptor_pdb": sha256_file(receptor_path),
        "ligands_sdf": sha256_file(ligand_path),
        "candidate_manifest_csv": sha256_file(candidate_manifest_path),
        "rejected_ligands_csv": sha256_file(rejected_path),
    }
    output_hashes["docking_config_json"] = sha256_file(output_dir / "docking_config.json")
    input_hashes = {
        "consensus_csv": sha256_file(consensus_path),
        "source_pdb": sha256_file(source_pdb),
    }
    if args.ranked_csv:
        input_hashes["ranked_csv"] = sha256_file(args.ranked_csv)
    if args.candidate_list:
        input_hashes["candidate_list_csv"] = sha256_file(args.candidate_list)
    write_json(
        output_dir / "input_run.json",
        {
            "schema_version": 2,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            **config,
            "input_sha256": input_hashes,
            "output_sha256": output_hashes,
        },
    )
    status_token = "GNINA_INPUTS_OK" if not rejected else "GNINA_INPUTS_PARTIAL"
    print(
        f"{status_token} pdb_id={pdb_id} consensus_id={consensus_id} "
        f"candidates={requested_count} ligands={len(accepted)} rejected={len(rejected)} box="
        f"{box['size_x']}x{box['size_y']}x{box['size_z']} purpose={args.purpose}"
    )
    print(output_dir)
    return 0 if accepted and not rejected else 2


if __name__ == "__main__":
    raise SystemExit(main())
