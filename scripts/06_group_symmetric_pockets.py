from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser(description="Group symmetry-equivalent P2Rank pockets by author residue numbers.")
    parser.add_argument("--pdb-id", required=True)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--jaccard-threshold", type=float, default=0.6)
    args = parser.parse_args()

    pdb_id = args.pdb_id.upper()
    source = ROOT / "runs" / "p2rank" / pdb_id / "output" / f"{pdb_id}.pdb_predictions.csv"
    with source.open(encoding="utf-8", newline="") as handle:
        raw_rows = [{key.strip(): value.strip() for key, value in row.items()} for row in csv.DictReader(handle)]
    pockets = []
    for row in raw_rows[: args.top]:
        residues = {token.split("_", 1)[1] for token in row["residue_ids"].split()}
        pockets.append(
            {
                "name": row["name"],
                "rank": int(row["rank"]),
                "probability": float(row["probability"]),
                "residues": residues,
            }
        )

    parent = list(range(len(pockets)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for left in range(len(pockets)):
        for right in range(left + 1, len(pockets)):
            a, b = pockets[left]["residues"], pockets[right]["residues"]
            jaccard = len(a & b) / len(a | b) if a or b else 0.0
            if jaccard >= args.jaccard_threshold:
                union(left, right)

    grouped: dict[int, list[dict[str, object]]] = {}
    for index, pocket in enumerate(pockets):
        grouped.setdefault(find(index), []).append(pocket)

    rows = []
    for group_index, members in enumerate(sorted(grouped.values(), key=lambda value: min(item["rank"] for item in value)), 1):
        residue_sets = [member["residues"] for member in members]
        union_residues = set().union(*residue_sets)
        core_residues = set.intersection(*residue_sets)
        rows.append(
            {
                "pdb_id": pdb_id,
                "group_id": f"G{group_index}",
                "member_count": len(members),
                "members": " ".join(member["name"] for member in sorted(members, key=lambda item: item["rank"])),
                "best_rank": min(member["rank"] for member in members),
                "max_probability": max(member["probability"] for member in members),
                "core_residue_count": len(core_residues),
                "union_residue_count": len(union_residues),
                "core_author_residue_numbers": " ".join(sorted(core_residues, key=lambda value: (int(''.join(filter(str.isdigit, value)) or 0), value))),
                "interpretation": "symmetry-equivalent candidate group; requires biological review",
            }
        )

    output_dir = ROOT / "results" / "pocket_groups"
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"{pdb_id}_symmetry_groups.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    json_path = output_dir / f"{pdb_id}_symmetry_groups.json"
    json_path.write_text(json.dumps(rows, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"{pdb_id}: grouped {len(pockets)} pockets into {len(rows)} residue-signature groups")
    print(f"Wrote {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
