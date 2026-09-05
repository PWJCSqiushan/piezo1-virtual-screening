from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.assessment import read_rows, write_rows
from piezo_vs.io_utils import sha256_file, write_json


MASTER_FIELDS = [
    "compound_id", "canonical_smiles", "source", "source_id", "source_url",
    "upstream_rank", "upstream_score", "docking_rank", "CNNscore",
    "minimizedAffinity", "swissadme_raw_file", "admetlab3_raw_file",
    "protox3_raw_file", "assessment_branch", "review_status", "reviewer", "review_note",
]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare a lossless Step 8 review table and tool inputs from ranked compounds."
    )
    parser.add_argument("--ranked-csv", type=Path, required=True)
    parser.add_argument("--docking-csv", type=Path)
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.top_n < 1:
        raise SystemExit("--top-n must be positive")
    if not args.ranked_csv.is_file():
        raise FileNotFoundError(args.ranked_csv)
    if args.docking_csv and not args.docking_csv.is_file():
        raise FileNotFoundError(args.docking_csv)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    ranked = read_rows(args.ranked_csv)[: args.top_n]
    if not ranked:
        raise ValueError("Ranked CSV is empty")
    docking_by_id = {}
    if args.docking_csv:
        docking_rows = read_rows(args.docking_csv)
        docking_by_id = {row.get("compound_id", ""): row for row in docking_rows}

    master = []
    seen = set()
    for row in ranked:
        compound_id = (row.get("compound_id") or "").strip()
        smiles = (row.get("canonical_smiles") or "").strip()
        if not compound_id or not smiles or compound_id in seen:
            raise ValueError(f"Missing/duplicate compound_id or SMILES: {compound_id!r}")
        seen.add(compound_id)
        dock = docking_by_id.get(compound_id, {})
        master.append({
            "compound_id": compound_id,
            "canonical_smiles": smiles,
            "source": row.get("source", ""),
            "source_id": row.get("source_id", ""),
            "source_url": row.get("source_url", ""),
            "upstream_rank": row.get("rank", ""),
            "upstream_score": row.get("drugclip_score", ""),
            "docking_rank": dock.get("docking_rank", ""),
            "CNNscore": dock.get("CNNscore", ""),
            "minimizedAffinity": dock.get("minimizedAffinity", ""),
            "swissadme_raw_file": "",
            "admetlab3_raw_file": "",
            "protox3_raw_file": "",
            "assessment_branch": "pending_tool_exports",
            "review_status": "pending",
            "reviewer": "",
            "review_note": "computational candidate; not efficacy or safety evidence",
        })

    master_path = args.output_dir / "compound_assessment_master.csv"
    tool_input_path = args.output_dir / "tool_input.csv"
    write_rows(master_path, master, MASTER_FIELDS)
    write_rows(
        tool_input_path,
        ({"compound_id": row["compound_id"], "canonical_smiles": row["canonical_smiles"]} for row in master),
        ["compound_id", "canonical_smiles"],
    )
    run = {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "prepared_waiting_for_external_predictions",
        "compound_count": len(master),
        "ranked_csv": str(args.ranked_csv.resolve()),
        "docking_csv": str(args.docking_csv.resolve()) if args.docking_csv else "",
        "scientific_status": "assessment_inputs_only_not_drug_safety_evidence",
        "input_sha256": {
            "ranked_csv": sha256_file(args.ranked_csv),
            "docking_csv": sha256_file(args.docking_csv) if args.docking_csv else "",
        },
        "output_sha256": {
            "master_csv": sha256_file(master_path),
            "tool_input_csv": sha256_file(tool_input_path),
        },
    }
    write_json(args.output_dir / "run.json", run)
    print(f"ASSESSMENT_INPUTS_OK compounds={len(master)} status={run['status']}")
    print(master_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
