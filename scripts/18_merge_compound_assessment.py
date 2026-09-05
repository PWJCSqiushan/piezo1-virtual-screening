from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.assessment import assessment_branch, index_tool_rows, normalized_key, prefixed_values, read_rows, write_rows
from piezo_vs.io_utils import sha256_file, write_json


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Merge raw Step 8 exports without discarding oral-heuristic failures."
    )
    parser.add_argument("--master", type=Path, required=True)
    parser.add_argument("--swissadme", type=Path)
    parser.add_argument("--admetlab3", type=Path)
    parser.add_argument("--protox3", type=Path)
    parser.add_argument(
        "--oral-pass-column",
        help="Exact SwissADME export column interpreted as yes/no oral-priority heuristic",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    inputs = {"swissadme": args.swissadme, "admetlab3": args.admetlab3, "protox3": args.protox3}
    if not args.master.is_file():
        raise FileNotFoundError(args.master)
    for path in inputs.values():
        if path and not path.is_file():
            raise FileNotFoundError(path)
    if args.oral_pass_column and not args.swissadme:
        raise SystemExit("--oral-pass-column requires --swissadme")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    master = read_rows(args.master)
    tool_indexes = {
        tool: index_tool_rows(read_rows(path), tool) if path else {}
        for tool, path in inputs.items()
    }
    merged = []
    for row in master:
        key = normalized_key(row)
        swiss = tool_indexes["swissadme"].get(key)
        branch, note = assessment_branch(swiss.get(args.oral_pass_column, "") if swiss and args.oral_pass_column else None)
        output = dict(row)
        output.update(prefixed_values(swiss, "swissadme"))
        output.update(prefixed_values(tool_indexes["admetlab3"].get(key), "admetlab3"))
        output.update(prefixed_values(tool_indexes["protox3"].get(key), "protox3"))
        output["swissadme_raw_file"] = str(args.swissadme.resolve()) if args.swissadme else ""
        output["admetlab3_raw_file"] = str(args.admetlab3.resolve()) if args.admetlab3 else ""
        output["protox3_raw_file"] = str(args.protox3.resolve()) if args.protox3 else ""
        output["assessment_branch"] = branch
        output["review_note"] = note
        merged.append(output)

    fields = list(dict.fromkeys(key for row in merged for key in row))
    assessed_path = args.output_dir / "compound_assessment_merged.csv"
    review_path = args.output_dir / "expert_review_queue.csv"
    write_rows(assessed_path, merged, fields)
    write_rows(review_path, (row for row in merged if row["assessment_branch"] != "oral_priority"), fields)
    write_json(args.output_dir / "run.json", {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "compound_count": len(merged),
        "branch_counts": {branch: sum(row["assessment_branch"] == branch for row in merged) for branch in sorted({row["assessment_branch"] for row in merged})},
        "oral_pass_column": args.oral_pass_column or "",
        "input_sha256": {"master": sha256_file(args.master), **{tool: sha256_file(path) if path else "" for tool, path in inputs.items()}},
        "output_sha256": {"merged_csv": sha256_file(assessed_path), "expert_review_queue_csv": sha256_file(review_path)},
        "scientific_status": "in_silico_triage_only_all_candidates_retained",
    })
    print(f"ASSESSMENT_MERGE_OK compounds={len(merged)} retained={len(merged)}")
    print(assessed_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
