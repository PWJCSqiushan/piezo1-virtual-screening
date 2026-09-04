from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.io_utils import read_json, sha256_file, write_json
from piezo_vs.project import ensure_project_dirs

USER_AGENT = "piezo1-vs/0.1 academic-competition"


def fetch(url: str, destination: Path, force: bool = False) -> dict[str, object]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not force:
        return {
            "url": url,
            "path": str(destination.relative_to(ROOT)),
            "bytes": destination.stat().st_size,
            "sha256": sha256_file(destination),
            "status": "cached",
        }

    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    partial = destination.with_suffix(destination.suffix + ".part")
    try:
        with urllib.request.urlopen(request, timeout=120) as response, partial.open("wb") as handle:
            while chunk := response.read(1024 * 1024):
                handle.write(chunk)
        partial.replace(destination)
        return {
            "url": url,
            "path": str(destination.relative_to(ROOT)),
            "bytes": destination.stat().st_size,
            "sha256": sha256_file(destination),
            "status": "downloaded",
        }
    except Exception:
        partial.unlink(missing_ok=True)
        raise


def fetch_json(url: str, destination: Path, force: bool = False) -> tuple[dict, dict[str, object]]:
    record = fetch(url, destination, force=force)
    return json.loads(destination.read_text(encoding="utf-8")), record


def main() -> int:
    parser = argparse.ArgumentParser(description="Download PIEZO1 structures and authoritative RCSB metadata.")
    parser.add_argument("--force", action="store_true", help="Redownload existing files.")
    args = parser.parse_args()

    ensure_project_dirs(ROOT)
    config = read_json(ROOT / "config" / "structures.json")
    records: list[dict[str, object]] = []
    errors: list[dict[str, str]] = []
    optional_warnings: list[dict[str, str]] = []

    for item in config["structures"]:
        pdb_id = item["pdb_id"].upper()
        lower = pdb_id.lower()
        raw_dir = ROOT / "data" / "raw" / "structures"
        api_dir = raw_dir / "api" / pdb_id

        targets = [
            (f"https://files.rcsb.org/download/{pdb_id}.cif", raw_dir / f"{pdb_id}.cif", True),
            (f"https://files.rcsb.org/download/{pdb_id}.pdb", raw_dir / f"{pdb_id}.pdb", True),
            (
                f"https://files.rcsb.org/pub/pdb/validation_reports/{lower[1:3]}/{lower}/{lower}_full_validation.pdf",
                raw_dir / f"{pdb_id}_full_validation.pdf",
                False,
            ),
        ]
        for url, destination, required in targets:
            try:
                records.append(fetch(url, destination, force=args.force))
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
                issue = {"pdb_id": pdb_id, "url": url, "error": str(exc)}
                (errors if required else optional_warnings).append(issue)

        entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        entry, entry_record = fetch_json(entry_url, api_dir / "entry.json", force=args.force)
        records.append(entry_record)

        entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
        for entity_id in entity_ids:
            url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
            try:
                _, record = fetch_json(url, api_dir / f"polymer_entity_{entity_id}.json", force=args.force)
                records.append(record)
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
                errors.append({"pdb_id": pdb_id, "url": url, "error": str(exc)})

    manifest = {
        "schema_version": 1,
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "records": records,
        "errors": errors,
        "optional_warnings": optional_warnings,
    }
    write_json(ROOT / "data" / "raw" / "structures" / "fetch_manifest.json", manifest)
    print(
        f"Downloaded/cached {len(records)} records; "
        f"errors={len(errors)}; optional_warnings={len(optional_warnings)}"
    )
    for error in errors:
        print(f"ERROR {error['pdb_id']} {error['url']}: {error['error']}", file=sys.stderr)
    for warning in optional_warnings:
        print(f"WARNING {warning['pdb_id']} {warning['url']}: {warning['error']}", file=sys.stderr)
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
