from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.io_utils import sha256_file, write_json


def read_ranks(path: Path) -> dict[str, int]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    ranks = {}
    for row in rows:
        compound_id = row.get("compound_id", "")
        rank = int(row.get("rank", "0"))
        if not compound_id or rank < 1 or compound_id in ranks:
            raise ValueError(f"Invalid or duplicate compound in {path}: {compound_id!r}")
        ranks[compound_id] = rank
    return ranks


def spearman_without_ties(first: dict[str, int], second: dict[str, int]) -> float:
    if first.keys() != second.keys():
        raise ValueError("DrugCLIP runs do not contain the same compound IDs")
    n = len(first)
    if n < 2:
        raise ValueError("At least two compounds are required")
    squared = sum((first[key] - second[key]) ** 2 for key in first)
    return 1.0 - (6.0 * squared) / (n * (n * n - 1))


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare two DrugCLIP rank lists for numerical stability.")
    parser.add_argument("--first", type=Path, required=True)
    parser.add_argument("--second", type=Path, required=True)
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.top_k < 1:
        raise SystemExit("--top-k must be positive")
    first = read_ranks(args.first)
    second = read_ranks(args.second)
    rho = spearman_without_ties(first, second)
    k = min(args.top_k, len(first))
    first_top = {key for key, rank in first.items() if rank <= k}
    second_top = {key for key, rank in second.items() if rank <= k}
    overlap = len(first_top & second_top)
    result = {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "compound_count": len(first),
        "spearman_rank_correlation": rho,
        "top_k": k,
        "top_k_overlap_count": overlap,
        "top_k_overlap_fraction": overlap / k,
        "identical_full_order": all(first[key] == second[key] for key in first),
        "inputs": {
            "first": str(args.first.resolve()),
            "first_sha256": sha256_file(args.first),
            "second": str(args.second.resolve()),
            "second_sha256": sha256_file(args.second),
        },
        "interpretation": "numerical_stability_check_only_not_biological_validation",
    }
    write_json(args.output, result)
    print(
        f"DRUGCLIP_STABILITY_OK compounds={len(first)} spearman={rho:.6f} "
        f"top{k}_overlap={overlap}/{k}"
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
