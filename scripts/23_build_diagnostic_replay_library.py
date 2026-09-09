from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.io_utils import sha256_file


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build a fixed diagnostic replay library from one preserved DrugCLIP export."
    )
    parser.add_argument("--task-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    task_dir = args.task_dir.resolve()
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output directory: {output_dir}")
    ranked_path = task_dir / "ranked_compounds.csv"
    manifest_path = task_dir / "molecule_manifest.csv"
    metadata_path = task_dir / "input_run.json"
    for path in (ranked_path, manifest_path, metadata_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    ranked = read_csv(ranked_path)
    manifest = read_csv(manifest_path)
    metadata = json.loads(metadata_path.read_text(encoding="utf-8-sig"))
    manifest_by_id = {row.get("compound_id", ""): row for row in manifest}
    if len(manifest_by_id) != len(manifest):
        raise ValueError("Source manifest has duplicate or empty compound IDs")
    output_rows: list[dict[str, str]] = []
    seen_smiles: set[str] = set()
    duplicate_smiles = 0
    for row in ranked:
        compound_id = row.get("compound_id", "")
        source = manifest_by_id.get(compound_id)
        if source is None or source.get("smiles") != row.get("smiles"):
            raise ValueError(f"Ranked row does not match manifest: {compound_id}")
        smiles = row.get("smiles", "")
        if smiles in seen_smiles:
            duplicate_smiles += 1
        seen_smiles.add(smiles)
        output_rows.append(
            {
                "compound_id": compound_id,
                "smiles": smiles,
                "source": "diagnostic_replay_from_teammate_export",
                "source_id": source.get("source_id", ""),
                "source_url": source.get("source_url", ""),
            }
        )
    if not output_rows:
        raise ValueError("Source task contains no ranked molecules")

    output_dir.mkdir(parents=True, exist_ok=True)
    library_path = output_dir / "diagnostic_library.csv"
    with library_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    provenance = {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "technical_validation_only",
        "formal_screening_eligible": False,
        "reason": (
            "The exported C006 rows are reused solely to force C006 and C007 through the exact "
            "same molecule universe. The original comp library identity, full size, license, "
            "and web score definition remain unverified."
        ),
        "source_task": {
            "path": str(task_dir),
            "pdb_id": metadata.get("pdb_id"),
            "consensus_id": metadata.get("consensus_id"),
            "ranked_sha256": sha256_file(ranked_path),
            "manifest_sha256": sha256_file(manifest_path),
            "input_run_sha256": sha256_file(metadata_path),
        },
        "row_count": len(output_rows),
        "unique_compound_id_count": len({row["compound_id"] for row in output_rows}),
        "duplicate_input_smiles_count": duplicate_smiles,
        "output": {
            "path": str(library_path),
            "sha256": sha256_file(library_path),
        },
        "retention_note": "No source rows were removed; LMDB preparation may separately reject invalid or duplicate canonical SMILES and records every rejection.",
    }
    (output_dir / "provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"DIAGNOSTIC_LIBRARY_OK rows={len(output_rows)} duplicate_smiles={duplicate_smiles}"
    )
    print(library_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())