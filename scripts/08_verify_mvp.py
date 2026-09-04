from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)
    print(f"PASS: {message}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify one independently analyzed PIEZO1 conformation.")
    parser.add_argument("--pdb-id", default="8YEZ")
    args = parser.parse_args()
    pdb_id = args.pdb_id.upper()

    manifest = csv_rows(ROOT / "results" / "structure_manifest.csv")
    rules = json.loads((ROOT / "config" / "consensus_rules.json").read_text(encoding="utf-8"))
    require(
        rules["active_pdb_ids"] == ["8YEZ", "8ZU3", "8YFC", "9VMX"],
        "active analysis scope contains exactly 8YEZ, 8ZU3, 8YFC and 9VMX",
    )
    require(rules["forbid_cross_conformation_consensus"] is True, "cross-conformation consensus is forbidden")
    require(rules["reference_only_pdb_ids"] == ["8ZU8"], "8ZU8 is reference-only")
    require(pdb_id in rules["active_pdb_ids"], f"{pdb_id} is an active independently analyzed conformation")
    require(any(row["pdb_id"] == pdb_id for row in manifest), f"{pdb_id} is present in the structure manifest")
    require(rules["tiers"]["T1"]["support_count"] == 2, "T1 means support from exactly two tools")
    require(rules["tiers"]["T2"]["support_count"] == 3, "T2 means support from all three tools")

    structure_paths = {}
    for suffix in ("pdb", "cif"):
        path = ROOT / "data" / "raw" / "structures" / f"{pdb_id}.{suffix}"
        require(path.exists() and path.stat().st_size > 1_000_000, f"{pdb_id}.{suffix} is present and non-empty")
        structure_paths[suffix] = path

    p2rank_path = ROOT / "runs" / "p2rank" / pdb_id / "output" / f"{pdb_id}.pdb_predictions.csv"
    p2rank_count = len(csv_rows(p2rank_path))
    require(p2rank_count > 0, f"{pdb_id} P2Rank output is present with {p2rank_count} pockets")

    dogsite_base = ROOT / "runs" / "dogsite3" / pdb_id
    latest = json.loads((dogsite_base / "latest_job.json").read_text(encoding="utf-8"))
    run_dir = dogsite_base / latest["job_id"]
    result = json.loads((run_dir / "result.json").read_text(encoding="utf-8"))
    require(result["status_code"] == 200, f"{pdb_id} DoGSite3 job completed")
    dogsite_count = len(result["residues"])
    require(dogsite_count > 0, f"{pdb_id} DoGSite3 output is present with {dogsite_count} pockets")
    downloads = json.loads((run_dir / "download_manifest.json").read_text(encoding="utf-8"))
    require(not [item for item in downloads if "error" in item], f"{pdb_id} DoGSite3 download manifest has no errors")
    require(all((ROOT / item["path"]).exists() for item in downloads), f"{pdb_id} DoGSite3 linked files exist locally")

    matches = csv_rows(ROOT / "results" / "pocket_matching" / f"{pdb_id}_p2rank_dogsite3_matches.csv")
    passed = sum(row["two_tool_candidate"] == "True" for row in matches)
    require(len(matches) == 20, f"{pdb_id} top-20 pairwise precheck is present; supported={passed}/20")

    print(f"PASS: approved analysis is restricted to one conformation at a time ({pdb_id})")
    fpocket_latest = ROOT / "runs" / "fpocket" / pdb_id / "latest_run.json"
    if fpocket_latest.exists():
        fpocket_run = json.loads(fpocket_latest.read_text(encoding="utf-8"))
        output_dir = ROOT / fpocket_run["output_dir"]
        pocket_files = list((output_dir / "pockets").glob("pocket*_atm.pdb"))
        require(fpocket_run["pdb_id"] == pdb_id, f"{pdb_id} fpocket metadata belongs to this structure")
        require(fpocket_run["returncode"] == 0, f"{pdb_id} fpocket process completed successfully")
        require(
            fpocket_run["input_sha256"] == sha256_file(structure_paths["pdb"]),
            f"{pdb_id} fpocket input SHA256 matches the configured PDB",
        )
        require(
            len(pocket_files) == fpocket_run["pocket_count"] > 0,
            f"{pdb_id} fpocket output is present with {len(pocket_files)} pockets",
        )
        consensus_path = ROOT / "results" / "consensus" / f"{pdb_id}_consensus_pockets.csv"
        consensus = csv_rows(consensus_path)
        require(consensus, f"{pdb_id} final three-tool consensus is present")
        require(
            all(row["support_count"] in {"2", "3"} and row["tier"] in {"T1", "T2"} for row in consensus),
            f"{pdb_id} final consensus contains only configured support_count 2/3 tiers",
        )
        tier1 = sum(row["support_count"] == "2" for row in consensus)
        tier2 = sum(row["support_count"] == "3" for row in consensus)
        print(f"INFO: tools completed 3/3; final consensus regions={len(consensus)}; T1={tier1}; T2={tier2}")
    else:
        print("INFO: tools completed 2/3; fpocket is pending")
        print("INFO: final T1/T2 classification has not been generated")
    print("SINGLE_CONFORMATION_STAGE_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
