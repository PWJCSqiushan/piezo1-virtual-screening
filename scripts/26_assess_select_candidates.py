from pathlib import Path
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from piezo_vs.candidate_workflow import build_snapshot


def main():
    parser = argparse.ArgumentParser(description="Evaluate an audited pool, import raw predictions, and select only approved candidates.")
    parser.add_argument("--batch-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--external-manifest", type=Path, action="append", default=[])
    parser.add_argument("--human-review", type=Path)
    args = parser.parse_args()
    status = build_snapshot(args.batch_dir.resolve(), args.output_dir.resolve(), args.external_manifest, args.human_review)
    print(json.dumps(status, indent=2))
    return 0 if status["selected_count"] == status["requested_count"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
