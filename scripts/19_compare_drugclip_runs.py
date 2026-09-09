from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.io_utils import sha256_file, write_json


def read_ranks(path: Path) -> tuple[dict[str, int], dict[str, float]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    ranks: dict[str, int] = {}
    scores: dict[str, float] = {}
    for row in rows:
        compound_id = row.get("compound_id", "")
        rank = int(row.get("rank", "0"))
        if not compound_id or rank < 1 or compound_id in ranks:
            raise ValueError(f"Invalid or duplicate compound in {path}: {compound_id!r}")
        ranks[compound_id] = rank
        scores[compound_id] = float(row.get("drugclip_score", "nan"))
    return ranks, scores


def spearman_without_ties(first: dict[str, int], second: dict[str, int]) -> float:
    if first.keys() != second.keys():
        raise ValueError("DrugCLIP runs do not contain the same compound IDs")
    n = len(first)
    if n < 2:
        raise ValueError("At least two compounds are required")
    squared = sum((first[key] - second[key]) ** 2 for key in first)
    return 1.0 - (6.0 * squared) / (n * (n * n - 1))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare two DrugCLIP rank lists only when they used the same compound universe."
    )
    parser.add_argument("--first", type=Path, required=True)
    parser.add_argument("--second", type=Path, required=True)
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--require-same-universe",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Return exit code 2 when the compound universes differ (default: true).",
    )
    args = parser.parse_args()
    if args.top_k < 1:
        raise SystemExit("--top-k must be positive")

    first, first_scores = read_ranks(args.first)
    second, second_scores = read_ranks(args.second)
    first_ids = set(first)
    second_ids = set(second)
    common = first_ids & second_ids
    union = first_ids | second_ids
    same_universe = first_ids == second_ids
    comparable = same_universe and len(first) >= 2
    k = min(args.top_k, len(first), len(second))
    first_top = {key for key, rank in first.items() if rank <= k}
    second_top = {key for key, rank in second.items() if rank <= k}
    overlap = len(first_top & second_top)
    rho = spearman_without_ties(first, second) if comparable else None

    reasons: list[str] = []
    if not same_universe:
        reasons.append("compound_universes_differ")
    if same_universe and len(first) < 2:
        reasons.append("fewer_than_two_compounds")
    result = {
        "schema_version": 2,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "comparison_status": "comparable" if comparable else "not_comparable",
        "not_comparable_reasons": reasons,
        "same_compound_universe": same_universe,
        "first_compound_count": len(first),
        "second_compound_count": len(second),
        "common_compound_count": len(common),
        "union_compound_count": len(union),
        "compound_universe_jaccard": len(common) / len(union) if union else 0.0,
        "spearman_rank_correlation": rho,
        "top_k": k,
        "top_k_overlap_count": overlap,
        "top_k_overlap_fraction": overlap / k if k else 0.0,
        "identical_full_order": comparable and all(first[key] == second[key] for key in first),
        "identical_scores": comparable and all(
            first_scores[key] == second_scores[key] for key in first
        ),
        "inputs": {
            "first": str(args.first.resolve()),
            "first_sha256": sha256_file(args.first),
            "second": str(args.second.resolve()),
            "second_sha256": sha256_file(args.second),
        },
        "interpretation": "same-library numerical comparison only; never biological validation",
    }
    write_json(args.output, result)
    if comparable:
        print(
            f"DRUGCLIP_COMPARISON_OK compounds={len(first)} spearman={rho:.6f} "
            f"top{k}_overlap={overlap}/{k}"
        )
    else:
        print(
            f"DRUGCLIP_NOT_COMPARABLE first={len(first)} second={len(second)} "
            f"common={len(common)} universe_jaccard={result['compound_universe_jaccard']:.6f}"
        )
    print(args.output)
    return 2 if args.require_same_universe and not comparable else 0


if __name__ == "__main__":
    raise SystemExit(main())