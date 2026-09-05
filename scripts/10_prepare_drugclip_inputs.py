from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.drugclip_io import (
    extract_pocket_atoms,
    parse_residue_tokens,
    read_compound_inputs,
    select_consensus_row,
    write_lmdb_records,
    write_pocket_pdb,
)
from piezo_vs.io_utils import read_json, sha256_file, write_json


def generate_molecule_record(smiles: str, seed: int) -> tuple[dict[str, object], str]:
    try:
        import numpy as np
        from rdkit import Chem, rdBase
        from rdkit.Chem import AllChem
    except ImportError as exc:
        raise RuntimeError(
            "Molecule preparation needs numpy and RDKit; run bootstrap_drugclip_wsl.ps1 first"
        ) from exc

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("RDKit rejected the SMILES")
    canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
    molecule = Chem.AddHs(molecule)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = seed
    parameters.numThreads = 1
    if AllChem.EmbedMolecule(molecule, parameters) != 0:
        raise ValueError("RDKit could not generate a 3D conformer")
    try:
        if AllChem.MMFFHasAllMoleculeParams(molecule):
            AllChem.MMFFOptimizeMolecule(molecule, maxIters=500)
        else:
            AllChem.UFFOptimizeMolecule(molecule, maxIters=500)
    except Exception:
        # A usable embedded conformer is preferable to discarding a molecule
        # solely because force-field optimization did not converge.
        pass
    molecule = Chem.RemoveHs(molecule)
    coordinates = np.asarray(molecule.GetConformer().GetPositions(), dtype=np.float32)
    record: dict[str, object] = {
        "atoms": [atom.GetSymbol() for atom in molecule.GetAtoms()],
        "coordinates": [coordinates],
        "smi": canonical_smiles,
    }
    return record, rdBase.rdkitVersion


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create one-pocket DrugCLIP LMDB inputs from a final same-structure consensus."
    )
    parser.add_argument("--pdb-id", required=True)
    parser.add_argument("--consensus-id", required=True)
    parser.add_argument("--compounds", type=Path, required=True, help="CSV/TSV/SMI with SMILES values")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument(
        "--residue-policy",
        choices=("two_plus", "all_supporting_tools"),
        default="two_plus",
        help=(
            "Consensus residues used for the pocket. two_plus follows Step 4 (default); "
            "all_supporting_tools is a stricter sensitivity-analysis option."
        ),
    )
    parser.add_argument(
        "--allow-non-piezo1-scope",
        action="store_true",
        help="Explicitly allow MDFIC-only or PIEZO1-MDFIC interface pockets after expert review.",
    )
    args = parser.parse_args()

    pdb_id = args.pdb_id.upper()
    consensus_id = args.consensus_id.upper()
    rules = read_json(ROOT / "config" / "consensus_rules.json")
    if pdb_id not in rules["active_pdb_ids"]:
        raise SystemExit(f"{pdb_id} is not an active independently analyzed structure")

    consensus_path = ROOT / "results" / "consensus" / f"{pdb_id}_consensus_pockets.csv"
    structure_path = ROOT / "data" / "raw" / "structures" / f"{pdb_id}.pdb"
    for required in (consensus_path, structure_path, args.compounds):
        if not required.is_file():
            raise FileNotFoundError(required)
    tier_by_support = {
        int(definition["support_count"]): tier for tier, definition in rules["tiers"].items()
    }
    row = select_consensus_row(
        consensus_path,
        pdb_id,
        consensus_id,
        tier_by_support=tier_by_support,
        required_tools=set(rules["required_tools"]),
    )
    if row.get("pocket_scope") != "piezo1" and not args.allow_non_piezo1_scope:
        raise ValueError(
            f"{pdb_id}/{consensus_id} has pocket_scope={row.get('pocket_scope')!r}. "
            "Only PIEZO1-only pockets enter DrugCLIP by default; use "
            "--allow-non-piezo1-scope only after structural/medical review."
        )
    residue_field = {
        "two_plus": "core_residues_2plus",
        "all_supporting_tools": "core_residues_all_supporting_tools",
    }[args.residue_policy]
    residues = parse_residue_tokens(row[residue_field])
    pocket_atoms, missing_residues = extract_pocket_atoms(structure_path, residues)
    if missing_residues:
        missing = " ".join(f"{chain}_{number}" for chain, number in sorted(missing_residues))
        raise ValueError(f"Consensus residues absent from {structure_path.name}: {missing}")
    if not pocket_atoms:
        raise ValueError("No protein atoms were extracted for the selected consensus pocket")

    started_at = datetime.now(timezone.utc)
    if args.output_dir:
        output_dir = args.output_dir.resolve()
    else:
        run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
        output_dir = ROOT / "runs" / "drugclip_inputs" / pdb_id / consensus_id / run_id
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    compounds = read_compound_inputs(args.compounds)
    if not compounds:
        raise ValueError(f"No compounds found in {args.compounds}")
    molecule_records: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []
    rejected_rows: list[dict[str, object]] = []
    seen_smiles: set[str] = set()
    rdkit_version = ""
    for compound in compounds:
        try:
            record, rdkit_version = generate_molecule_record(compound.smiles, args.seed)
            canonical_smiles = str(record["smi"])
            if canonical_smiles in seen_smiles:
                rejected_rows.append(
                    {
                        "compound_id": compound.compound_id,
                        "input_smiles": compound.smiles,
                        "source_row": compound.source_row,
                        "reason": "duplicate_canonical_smiles",
                    }
                )
                continue
            seen_smiles.add(canonical_smiles)
            molecule_records.append(record)
            manifest_rows.append(
                {
                    "lmdb_index": len(molecule_records) - 1,
                    "compound_id": compound.compound_id,
                    "input_smiles": compound.smiles,
                    "canonical_smiles": canonical_smiles,
                    "source_row": compound.source_row,
                    "source": compound.source,
                    "source_id": compound.source_id,
                    "source_url": compound.source_url,
                    "atom_count": len(record["atoms"]),
                }
            )
        except ValueError as exc:
            rejected_rows.append(
                {
                    "compound_id": compound.compound_id,
                    "input_smiles": compound.smiles,
                    "source_row": compound.source_row,
                    "reason": str(exc),
                }
            )
    if not molecule_records:
        raise ValueError("RDKit could not prepare any valid compounds")

    pocket_pdb = output_dir / "standardized_pocket.pdb"
    write_pocket_pdb(
        pocket_pdb,
        pocket_atoms,
        pdb_id=pdb_id,
        consensus_id=consensus_id,
        support_count=int(row["support_count"]),
    )

    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("Pocket preparation needs numpy") from exc
    pocket_record = {
        "pocket": f"{pdb_id}:{consensus_id}",
        "pocket_atoms": [atom.atom_name for atom in pocket_atoms],
        "pocket_coordinates": [np.asarray(atom.coordinate, dtype=np.float32) for atom in pocket_atoms],
    }
    mol_lmdb = output_dir / "mols.lmdb"
    pocket_lmdb = output_dir / "pocket.lmdb"
    write_lmdb_records(mol_lmdb, molecule_records)
    write_lmdb_records(pocket_lmdb, [pocket_record])

    manifest_path = output_dir / "molecule_manifest.csv"
    with manifest_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]))
        writer.writeheader()
        writer.writerows(manifest_rows)
    rejected_path = output_dir / "rejected_compounds.csv"
    with rejected_path.open("w", encoding="utf-8-sig", newline="") as handle:
        fields = ["compound_id", "input_smiles", "source_row", "reason"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rejected_rows)

    metadata = {
        "schema_version": 1,
        "created_at_utc": started_at.isoformat(),
        "pdb_id": pdb_id,
        "consensus_id": consensus_id,
        "tier": row["tier"],
        "support_count": int(row["support_count"]),
        "tools": row["tools"],
        "pocket_scope": row["pocket_scope"],
        "residue_source_field": residue_field,
        "residue_policy": args.residue_policy,
        "residue_count": len(residues),
        "pocket_atom_count": len(pocket_atoms),
        "model_max_pocket_atoms": 256,
        "model_atom_crop_expected": len(pocket_atoms) > 256,
        "input_compound_count": len(compounds),
        "accepted_compound_count": len(molecule_records),
        "rejected_compound_count": len(rejected_rows),
        "seed": args.seed,
        "rdkit_version": rdkit_version,
        "source_files": {
            "structure": str(structure_path.relative_to(ROOT)),
            "structure_sha256": sha256_file(structure_path),
            "consensus": str(consensus_path.relative_to(ROOT)),
            "consensus_sha256": sha256_file(consensus_path),
            "compounds": str(args.compounds.resolve()),
            "compounds_sha256": sha256_file(args.compounds),
        },
        "outputs": {
            "mols_lmdb_sha256": sha256_file(mol_lmdb),
            "pocket_lmdb_sha256": sha256_file(pocket_lmdb),
            "molecule_manifest_sha256": sha256_file(manifest_path),
            "standardized_pocket_pdb": str(pocket_pdb.relative_to(ROOT)),
            "standardized_pocket_pdb_sha256": sha256_file(pocket_pdb),
        },
    }
    write_json(output_dir / "input_run.json", metadata)
    print(
        f"DRUGCLIP_INPUTS_OK pdb_id={pdb_id} consensus_id={consensus_id} tier={row['tier']} "
        f"molecules={len(molecule_records)} rejected={len(rejected_rows)} "
        f"pocket_atoms={len(pocket_atoms)}"
    )
    print(output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
