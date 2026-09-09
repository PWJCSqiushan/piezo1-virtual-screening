from __future__ import annotations

import argparse
import csv
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
USER_AGENT = "piezo1-vs/0.2 academic-competition structure-verification"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().lower()


def atom_record_hash(path: Path) -> tuple[str, int]:
    lines: list[str] = []
    with path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")):
                lines.append(line.rstrip("\r\n"))
    payload = ("\n".join(lines) + "\n").encode("ascii", errors="replace")
    return hashlib.sha256(payload).hexdigest().lower(), len(lines)


def read_manifest(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return {row["pdb_id"].upper(): row for row in csv.DictReader(handle)}


def download(url: str, destination: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=120) as response:
        data = response.read()
    if not data.startswith(b"HEADER"):
        raise ValueError(f"RCSB response for {url} is not a legacy PDB file")
    destination.write_bytes(data)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fresh-download RCSB PDB files into an isolated cache and compare atom hashes."
    )
    parser.add_argument("--pdb-ids", nargs="+", default=["8YEZ", "8ZU3", "8YFC", "9VMX"])
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = ROOT / "results" / "structure_manifest.csv"
    manifest = read_manifest(manifest_path)
    rows: list[dict[str, object]] = []

    for raw_id in args.pdb_ids:
        pdb_id = raw_id.upper()
        if pdb_id not in manifest:
            raise ValueError(f"{pdb_id} is absent from structure_manifest.csv")
        local_path = ROOT / "data" / "raw" / "structures" / f"{pdb_id}.pdb"
        if not local_path.is_file():
            raise FileNotFoundError(local_path)
        fresh_path = output_dir / f"{pdb_id}.pdb"
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        download(url, fresh_path)
        local_atom_hash, local_atom_count = atom_record_hash(local_path)
        fresh_atom_hash, fresh_atom_count = atom_record_hash(fresh_path)
        expected = manifest[pdb_id]
        rows.append(
            {
                "pdb_id": pdb_id,
                "rcsb_url": url,
                "local_path": str(local_path),
                "fresh_path": str(fresh_path),
                "local_file_sha256": sha256(local_path),
                "fresh_file_sha256": sha256(fresh_path),
                "whole_file_exact_match": sha256(local_path) == sha256(fresh_path),
                "local_atom_sha256": local_atom_hash,
                "fresh_atom_sha256": fresh_atom_hash,
                "atom_records_exact_match": local_atom_hash == fresh_atom_hash,
                "local_atom_count": local_atom_count,
                "fresh_atom_count": fresh_atom_count,
                "manifest_atom_sha256": expected["coordinate_atom_sha256"],
                "manifest_atom_hash_matches_local": expected["coordinate_atom_sha256"] == local_atom_hash,
                "manifest_pdb_sha256": expected["pdb_sha256"],
                "manifest_file_hash_matches_local": expected["pdb_sha256"] == sha256(local_path),
                "mutation_label": expected["mutation_label"],
                "mutation_coordinate_status": expected["mutation_coordinate_status"],
            }
        )

    fields = list(rows[0])
    csv_path = output_dir / "rcsb_structure_verification.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    by_atom_hash: dict[str, list[str]] = {}
    for row in rows:
        by_atom_hash.setdefault(str(row["fresh_atom_sha256"]), []).append(str(row["pdb_id"]))
    identical_groups = [group for group in by_atom_hash.values() if len(group) > 1]
    result = {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": "fresh HTTPS download from files.rcsb.org",
        "all_fresh_match_local_atom_records": all(
            bool(row["atom_records_exact_match"]) for row in rows
        ),
        "all_local_files_match_manifest": all(
            bool(row["manifest_file_hash_matches_local"]) for row in rows
        ),
        "fresh_identical_atom_groups": identical_groups,
        "interpretation": (
            "If fresh RCSB files match local files, identical modeled atom records are present "
            "in the official deposits rather than introduced by this repository. Mutation "
            "biology still requires expert interpretation, especially when the mutation site "
            "is absent from modeled coordinates."
        ),
        "rows": rows,
    }
    (output_dir / "rcsb_structure_verification.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"RCSB_STRUCTURE_VERIFY_OK structures={len(rows)} "
        f"fresh_match_local={result['all_fresh_match_local_atom_records']} "
        f"identical_groups={identical_groups}"
    )
    print(output_dir)
    return 0 if result["all_fresh_match_local_atom_records"] else 2


if __name__ == "__main__":
    raise SystemExit(main())