from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def clean_row(row: dict[str, str]) -> dict[str, str]:
    return {key.strip(): value.strip() for key, value in row.items()}


def parse_p2rank(path: Path) -> list[dict[str, object]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [clean_row(row) for row in csv.DictReader(handle)]
    pockets = []
    for row in rows:
        residues = set()
        for token in row["residue_ids"].split():
            chain, residue_number = token.split("_", 1)
            residues.add((chain, residue_number))
        pockets.append(
            {
                "name": row["name"],
                "rank": int(row["rank"]),
                "score": float(row["score"]),
                "probability": float(row["probability"]),
                "center": tuple(float(row[key]) for key in ("center_x", "center_y", "center_z")),
                "residues": residues,
            }
        )
    return pockets


def parse_dogsite_residues(path: Path) -> tuple[set[tuple[str, str]], tuple[float, float, float]]:
    residues: set[tuple[str, str]] = set()
    coordinates: list[tuple[float, float, float]] = []
    with path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            chain = line[21].strip() or "_"
            residue_number = line[22:27].strip()
            residues.add((chain, residue_number))
            coordinates.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not coordinates:
        raise ValueError(f"No atom coordinates in {path}")
    count = len(coordinates)
    center = tuple(sum(point[index] for point in coordinates) / count for index in range(3))
    return residues, center


def parse_dogsite_table(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {row["name"]: row for row in csv.DictReader(handle, delimiter="\t")}


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((left - right) ** 2 for left, right in zip(a, b)))


def latest_dogsite_run(pdb_id: str) -> Path:
    base = ROOT / "runs" / "dogsite3" / pdb_id
    latest = json.loads((base / "latest_job.json").read_text(encoding="utf-8"))
    return base / latest["job_id"]


def main() -> int:
    parser = argparse.ArgumentParser(description="Match P2Rank and DoGSite3 pockets by coordinates and residues.")
    parser.add_argument("--pdb-id", required=True)
    parser.add_argument("--top-p2rank", type=int, default=20)
    parser.add_argument("--max-center-distance", type=float, default=12.0)
    parser.add_argument("--min-shared-residues", type=int, default=3)
    args = parser.parse_args()

    pdb_id = args.pdb_id.upper()
    p2rank_path = ROOT / "runs" / "p2rank" / pdb_id / "output" / f"{pdb_id}.pdb_predictions.csv"
    dogsite_run = latest_dogsite_run(pdb_id)
    downloads = dogsite_run / "downloads"
    dogsite_table_path = downloads / f"{pdb_id.lower()}_desc.txt"

    p2rank = parse_p2rank(p2rank_path)[: args.top_p2rank]
    descriptors = parse_dogsite_table(dogsite_table_path)
    dogsite = []
    for residue_path in sorted(downloads.glob(f"{pdb_id.lower()}_P_*_res.pdb")):
        name = residue_path.name.removeprefix(f"{pdb_id.lower()}_").removesuffix("_res.pdb")
        residues, center = parse_dogsite_residues(residue_path)
        dogsite.append({"name": name, "residues": residues, "center": center, "descriptor": descriptors[name]})

    matches = []
    for p2 in p2rank:
        candidates = []
        for dog in dogsite:
            shared = p2["residues"] & dog["residues"]
            union = p2["residues"] | dog["residues"]
            overlap_denominator = min(len(p2["residues"]), len(dog["residues"]))
            center_distance = distance(p2["center"], dog["center"])
            jaccard = len(shared) / len(union) if union else 0.0
            overlap = len(shared) / overlap_denominator if overlap_denominator else 0.0
            candidates.append((len(shared), overlap, jaccard, -center_distance, dog, center_distance))
        nearby = [candidate for candidate in candidates if candidate[5] <= args.max_center_distance]
        # First respect spatial proximity; residue numbering alone can spuriously match
        # distant regions in a very large, repeated or symmetric protein.
        shared_count, overlap, jaccard, _, dog, center_distance = max(nearby or candidates)
        descriptor = dog["descriptor"]
        matches.append(
            {
                "pdb_id": pdb_id,
                "p2rank_pocket": p2["name"],
                "p2rank_rank": p2["rank"],
                "p2rank_score": p2["score"],
                "p2rank_probability": p2["probability"],
                "p2rank_residue_count": len(p2["residues"]),
                "p2rank_chains": ",".join(sorted({item[0] for item in p2["residues"]})),
                "dogsite_pocket": dog["name"],
                "dogsite_residue_count": len(dog["residues"]),
                "dogsite_chains": ",".join(sorted({item[0] for item in dog["residues"]})),
                "dogsite_volume": descriptor.get("volume", ""),
                "dogsite_enclosure": descriptor.get("enclosure", ""),
                "dogsite_depth": descriptor.get("depth", ""),
                "dogsite_hydrophobicity": descriptor.get("hydrophobicity", ""),
                "center_distance_angstrom": round(center_distance, 3),
                "shared_residue_count": shared_count,
                "residue_jaccard": round(jaccard, 4),
                "residue_overlap_coefficient": round(overlap, 4),
                "two_tool_candidate": center_distance <= args.max_center_distance
                and shared_count >= args.min_shared_residues,
                "shared_residues": " ".join(f"{chain}_{number}" for chain, number in sorted(p2["residues"] & dog["residues"])),
            }
        )

    output_dir = ROOT / "results" / "pocket_matching"
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"{pdb_id}_p2rank_dogsite3_matches.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(matches[0]))
        writer.writeheader()
        writer.writerows(matches)
    json_path = output_dir / f"{pdb_id}_p2rank_dogsite3_matches.json"
    json_path.write_text(json.dumps(matches, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    passed = sum(bool(row["two_tool_candidate"]) for row in matches)
    print(f"{pdb_id}: matched {len(matches)} top P2Rank pockets; two-tool candidates={passed}")
    print(f"Wrote {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
