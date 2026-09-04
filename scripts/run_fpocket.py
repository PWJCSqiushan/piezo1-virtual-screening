from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.io_utils import read_json, sha256_file, write_json


def main() -> int:
    parser = argparse.ArgumentParser(description="Run fpocket for one configured PIEZO1 structure.")
    parser.add_argument("--pdb-id", required=True)
    parser.add_argument("--binary", default=os.environ.get("FPOCKET_BIN", "fpocket"))
    args = parser.parse_args()

    pdb_id = args.pdb_id.upper()
    rules = read_json(ROOT / "config" / "consensus_rules.json")
    if pdb_id not in rules["active_pdb_ids"]:
        raise SystemExit(f"{pdb_id} is not an active independently analyzed structure")
    binary = shutil.which(args.binary) if not Path(args.binary).is_file() else str(Path(args.binary).resolve())
    if not binary:
        raise SystemExit(
            "fpocket executable not found. Install fpocket in Linux/WSL/HPC, add it to PATH, "
            "or pass --binary /absolute/path/to/fpocket."
        )

    source = ROOT / "data" / "raw" / "structures" / f"{pdb_id}.pdb"
    if not source.exists():
        raise FileNotFoundError(f"Missing {source}; run scripts/01_fetch_structures.py first")

    started_at = datetime.now(timezone.utc)
    run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    run_dir = ROOT / "runs" / "fpocket" / pdb_id / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    run_input = run_dir / f"{pdb_id}.pdb"
    shutil.copy2(source, run_input)
    command = [str(binary), "-f", str(run_input)]
    process = subprocess.run(command, cwd=run_dir, text=True, capture_output=True, check=False)
    (run_dir / "stdout.log").write_text(process.stdout, encoding="utf-8")
    (run_dir / "stderr.log").write_text(process.stderr, encoding="utf-8")
    output_dir = run_dir / f"{pdb_id}_out"
    version_probe = subprocess.run([str(binary), "-h"], text=True, capture_output=True, check=False)
    version_text = (version_probe.stdout or version_probe.stderr).splitlines()
    metadata = {
        "schema_version": 1,
        "pdb_id": pdb_id,
        "run_id": run_id,
        "started_at_utc": started_at.isoformat(),
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "binary": str(binary),
        "version_banner": version_text[0].strip() if version_text else "",
        "command": command,
        "returncode": process.returncode,
        "input_sha256": sha256_file(source),
        "run_dir": str(run_dir.relative_to(ROOT)),
        "output_dir": str(output_dir.relative_to(ROOT)),
    }
    write_json(run_dir / "run.json", metadata)
    if process.returncode:
        print(process.stderr, file=sys.stderr)
        return process.returncode
    if not list((output_dir / "pockets").glob("pocket*_atm.pdb")):
        raise RuntimeError(f"fpocket returned success but no pocket*_atm.pdb files were found in {output_dir}")
    metadata["pocket_count"] = len(list((output_dir / "pockets").glob("pocket*_atm.pdb")))
    write_json(run_dir / "run.json", metadata)
    write_json(ROOT / "runs" / "fpocket" / pdb_id / "latest_run.json", metadata)
    print(f"fpocket output: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
