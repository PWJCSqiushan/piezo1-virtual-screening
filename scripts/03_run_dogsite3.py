from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.io_utils import sha256_file, write_json

API = "https://proteins.plus/api/dogsite3_rest"
USER_AGENT = "piezo1-vs/0.1 academic-competition"


def request_json(url: str, payload: dict[str, Any] | None = None) -> dict[str, Any]:
    data = json.dumps(payload).encode("utf-8") if payload is not None else None
    request = urllib.request.Request(
        url,
        data=data,
        headers={"Accept": "application/json", "Content-Type": "application/json", "User-Agent": USER_AGENT},
        method="POST" if payload is not None else "GET",
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def iter_urls(value: Any, label: str = "result"):
    if isinstance(value, str) and value.startswith(("http://", "https://")):
        yield label, value
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from iter_urls(item, f"{label}_{index + 1}")
    elif isinstance(value, dict):
        for key, item in value.items():
            yield from iter_urls(item, str(key))


def safe_filename(label: str, url: str) -> str:
    path_name = Path(urllib.parse.urlparse(url).path).name
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", path_name or label)
    return cleaned if cleaned else f"{label}.dat"


def download(url: str, destination: Path, retries: int = 3) -> dict[str, Any]:
    if destination.exists() and destination.stat().st_size > 0:
        return {
            "url": url,
            "path": str(destination.relative_to(ROOT)),
            "bytes": destination.stat().st_size,
            "sha256": sha256_file(destination),
            "status": "cached",
        }
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=180) as response, partial.open("wb") as handle:
                handle.write(response.read())
            partial.replace(destination)
            break
        except Exception:
            partial.unlink(missing_ok=True)
            if attempt == retries:
                raise
            time.sleep(attempt * 2)
    return {
        "url": url,
        "path": str(destination.relative_to(ROOT)),
        "bytes": destination.stat().st_size,
        "sha256": sha256_file(destination),
        "status": "downloaded",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Submit and collect a DoGSite3 pocket-prediction job.")
    parser.add_argument("--pdb-id", required=True)
    parser.add_argument("--chain", default="", help="Optional author chain; empty means all chains.")
    parser.add_argument("--poll-seconds", type=int, default=10)
    parser.add_argument("--timeout-seconds", type=int, default=1800)
    parser.add_argument("--resume-location", default="")
    args = parser.parse_args()

    pdb_id = args.pdb_id.upper()
    base_dir = ROOT / "runs" / "dogsite3" / pdb_id
    base_dir.mkdir(parents=True, exist_ok=True)

    if args.resume_location:
        location = args.resume_location
        submit_response = {"location": location, "status_code": 202, "resumed": True}
    else:
        payload = {
            "dogsite3": {
                "pdbCode": pdb_id,
                "analysisDetail": "0",
                "bindingSitePredictionGranularity": "1",
                "ligand": "",
                "chain": args.chain,
                "ligandBias": "0",
            }
        }
        submit_response = request_json(API, payload)
        location = submit_response["location"]

    job_id = location.rstrip("/").split("/")[-1]
    run_dir = base_dir / job_id
    run_dir.mkdir(parents=True, exist_ok=True)
    write_json(run_dir / "submit_response.json", submit_response)
    write_json(
        base_dir / "latest_job.json",
        {"pdb_id": pdb_id, "job_id": job_id, "location": location, "submitted_at_utc": datetime.now(timezone.utc).isoformat()},
    )
    print(f"DoGSite3 job {job_id}: {location}", flush=True)

    deadline = time.monotonic() + args.timeout_seconds
    response: dict[str, Any]
    while True:
        response = request_json(location)
        write_json(run_dir / "latest_response.json", response)
        status = int(response.get("status_code", 0))
        print(f"status={status} message={response.get('message', '')}", flush=True)
        if status == 200:
            break
        if status not in (201, 202):
            raise RuntimeError(f"Unexpected DoGSite3 response: {response}")
        if time.monotonic() >= deadline:
            raise TimeoutError(f"DoGSite3 job did not finish within {args.timeout_seconds}s; resume with {location}")
        time.sleep(args.poll_seconds)

    write_json(run_dir / "result.json", response)
    downloads = []
    used_names: set[str] = set()
    for label, url in iter_urls(response):
        if url.rstrip("/") == location.rstrip("/"):
            continue
        filename = safe_filename(label, url)
        if filename in used_names:
            filename = f"{label}_{filename}"
        used_names.add(filename)
        try:
            downloads.append(download(url, run_dir / "downloads" / filename))
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            downloads.append({"url": url, "error": str(exc)})
    write_json(run_dir / "download_manifest.json", downloads)
    failed = [item for item in downloads if "error" in item]
    print(
        f"Completed {pdb_id}; downloaded {sum('path' in item for item in downloads)} result files; "
        f"failed={len(failed)}"
    )
    if failed:
        print(f"ERROR: incomplete download; rerun with --resume-location {location}", file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
