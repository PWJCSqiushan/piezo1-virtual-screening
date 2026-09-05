from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.drugclip_io import parse_ranked_compounds
from piezo_vs.io_utils import sha256_file, write_json


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert DrugCLIP ranked text to an attributed CSV.")
    parser.add_argument("--ranked", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pdb-id", required=True)
    parser.add_argument("--consensus-id", required=True)
    parser.add_argument("--tier", choices=("T1", "T2"), required=True)
    args = parser.parse_args()

    with args.manifest.open(encoding="utf-8-sig", newline="") as handle:
        manifest_rows = list(csv.DictReader(handle))
    by_smiles = {row["canonical_smiles"]: row for row in manifest_rows}
    ranked = parse_ranked_compounds(args.ranked)
    missing = sorted({smiles for smiles, _ in ranked if smiles not in by_smiles})
    if missing:
        raise ValueError(f"{len(missing)} ranked SMILES are absent from the molecule manifest")
    output_rows = []
    for rank, (smiles, score) in enumerate(ranked, 1):
        source = by_smiles[smiles]
        output_rows.append(
            {
                "pdb_id": args.pdb_id.upper(),
                "consensus_id": args.consensus_id.upper(),
                "tier": args.tier,
                "rank": rank,
                "compound_id": source["compound_id"],
                "canonical_smiles": smiles,
                "drugclip_score": score,
                "interpretation": "virtual_screening_rank_only_not_experimental_evidence",
            }
        )
    if not output_rows:
        raise ValueError("DrugCLIP ranked output is empty")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite {args.output}")
    with args.output.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    write_json(
        args.output.with_suffix(args.output.suffix + ".run.json"),
        {
            "schema_version": 1,
            "pdb_id": args.pdb_id.upper(),
            "consensus_id": args.consensus_id.upper(),
            "tier": args.tier,
            "ranked_input_sha256": sha256_file(args.ranked),
            "manifest_input_sha256": sha256_file(args.manifest),
            "output_sha256": sha256_file(args.output),
            "row_count": len(output_rows),
            "scientific_status": "computational_prediction_only",
        },
    )
    print(f"DRUGCLIP_RESULTS_OK rows={len(output_rows)} output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

