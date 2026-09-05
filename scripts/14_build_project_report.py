from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.reporting import build_snapshot, write_report_bundle


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a self-contained project evidence report.")
    parser.add_argument("--ranked-csv", type=Path)
    parser.add_argument("--input-run", type=Path)
    parser.add_argument("--docking-run", type=Path)
    parser.add_argument("--assessment-run", type=Path)
    parser.add_argument("--stability-json", type=Path)
    parser.add_argument(
        "--purpose",
        choices=("technical_validation", "formal_screening"),
        default="technical_validation",
    )
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    if bool(args.ranked_csv) != bool(args.input_run):
        raise SystemExit("--ranked-csv and --input-run must be supplied together")
    for path in (args.ranked_csv, args.input_run, args.docking_run, args.assessment_run, args.stability_json):
        if path and not path.is_file():
            raise FileNotFoundError(path)
    output_dir = args.output_dir
    if output_dir is None:
        run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        output_dir = ROOT / "results" / "evidence_reports" / run_id
    snapshot = build_snapshot(
        ROOT,
        ranked_csv=args.ranked_csv,
        input_run=args.input_run,
        docking_run=args.docking_run,
        assessment_run=args.assessment_run,
        stability_json=args.stability_json,
        purpose=args.purpose,
    )
    outputs = write_report_bundle(output_dir.resolve(), snapshot)
    print(
        f"PROJECT_REPORT_OK structures={len(snapshot['structures'])} "
        f"consensus_regions={snapshot['total_consensus_regions']} "
        f"drugclip_rows={snapshot['drugclip']['ranked_count']}"
    )
    print(outputs["html"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
