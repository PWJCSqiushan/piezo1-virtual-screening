from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.docking import docking_box_from_atoms, read_ranked_compounds, write_receptor_pdb
from piezo_vs.drugclip_io import extract_pocket_atoms, parse_residue_tokens, select_consensus_row
from piezo_vs.io_utils import read_json, sha256_file, write_json


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


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare a traceable GNINA receptor, ligand SDF and pocket-centered search box."
    )
    parser.add_argument("--pdb-id", required=True)
    parser.add_argument("--consensus-id", required=True)
    parser.add_argument("--ranked-csv", type=Path, required=True)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--allow-non-piezo1-scope", action="store_true")
    parser.add_argument("--allow-large-box-clipping", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.top_n <= 100:
        raise SystemExit("top-n must be between 1 and 100")

    pdb_id = args.pdb_id.upper()
    consensus_id = args.consensus_id.upper()
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    structure_config = read_json(ROOT / "config" / "structures.json")
    structure = next(item for item in structure_config["structures"] if item["pdb_id"] == pdb_id)
    consensus_path = ROOT / "results" / "consensus" / f"{pdb_id}_consensus_pockets.csv"
    source_pdb = ROOT / "data" / "raw" / "structures" / f"{pdb_id}.pdb"
    for path in (args.ranked_csv, consensus_path, source_pdb):
        if not path.is_file():
            raise FileNotFoundError(path)
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

    residue_field = "core_residues_2plus"
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

    receptor_chains = set(structure["auth_piezo1_chains"])
    if scope in {"piezo1_mdfic_interface", "mdfic"}:
        receptor_chains.update(structure["auth_mdfic_chains"])
    receptor_path = output_dir / "receptor.pdb"
    receptor_atom_count = write_receptor_pdb(source_pdb, receptor_path, receptor_chains)

    rows = read_ranked_compounds(args.ranked_csv, pdb_id, consensus_id, args.top_n)
    try:
        from rdkit import Chem, rdBase
    except ImportError as exc:
        raise RuntimeError("GNINA input preparation requires RDKit") from exc
    ligand_path = output_dir / "ligands.sdf"
    writer = Chem.SDWriter(str(ligand_path))
    accepted: list[dict[str, str]] = []
    rejected: list[dict[str, str]] = []
    try:
        for index, row in enumerate(rows):
            try:
                molecule = prepare_ligand(row["canonical_smiles"], args.seed + index)
                molecule.SetProp("_Name", row["compound_id"])
                for key in (
                    "compound_id",
                    "canonical_smiles",
                    "source",
                    "source_id",
                    "source_url",
                    "rank",
                    "drugclip_score",
                ):
                    molecule.SetProp(key, row.get(key, ""))
                writer.write(molecule)
                accepted.append(row)
            except (ValueError, RuntimeError) as exc:
                rejected.append(
                    {
                        "compound_id": row.get("compound_id", ""),
                        "canonical_smiles": row.get("canonical_smiles", ""),
                        "reason": str(exc),
                    }
                )
    finally:
        writer.close()
    if not accepted:
        raise ValueError("No ranked compounds could be prepared for GNINA")

    rejected_path = output_dir / "rejected_ligands.csv"
    with rejected_path.open("w", encoding="utf-8-sig", newline="") as handle:
        fields = ["compound_id", "canonical_smiles", "reason"]
        csv_writer = csv.DictWriter(handle, fieldnames=fields)
        csv_writer.writeheader()
        csv_writer.writerows(rejected)

    config = {
        "schema_version": 1,
        "pdb_id": pdb_id,
        "consensus_id": consensus_id,
        "tier": consensus["tier"],
        "support_count": int(consensus["support_count"]),
        "pocket_scope": scope,
        "residue_source_field": residue_field,
        "box": box,
        "receptor": "receptor.pdb",
        "ligands": "ligands.sdf",
        "receptor_chains": sorted(receptor_chains),
        "receptor_atom_count": receptor_atom_count,
        "requested_top_n": args.top_n,
        "prepared_ligand_count": len(accepted),
        "rejected_ligand_count": len(rejected),
        "seed": args.seed,
        "rdkit_version": rdBase.rdkitVersion,
        "scientific_status": "computational_docking_input_not_binding_evidence",
    }
    write_json(output_dir / "docking_config.json", config)
    write_json(
        output_dir / "input_run.json",
        {
            "schema_version": 1,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            **config,
            "input_sha256": {
                "ranked_csv": sha256_file(args.ranked_csv),
                "consensus_csv": sha256_file(consensus_path),
                "source_pdb": sha256_file(source_pdb),
            },
            "output_sha256": {
                "receptor_pdb": sha256_file(receptor_path),
                "ligands_sdf": sha256_file(ligand_path),
                "rejected_ligands_csv": sha256_file(rejected_path),
            },
        },
    )
    print(
        f"GNINA_INPUTS_OK pdb_id={pdb_id} consensus_id={consensus_id} "
        f"ligands={len(accepted)} rejected={len(rejected)} box="
        f"{box['size_x']}x{box['size_y']}x{box['size_z']}"
    )
    print(output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
