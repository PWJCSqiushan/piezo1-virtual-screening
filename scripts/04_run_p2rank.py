from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.io_utils import read_json, sha256_file, write_json


def discover_java_home() -> str:
    configured = os.environ.get("JAVA_HOME")
    if configured:
        java_bin = Path(configured) / "bin" / ("java.exe" if os.name == "nt" else "java")
        if java_bin.exists():
            return configured
    probe = subprocess.run(
        ["java", "-XshowSettings:properties", "-version"],
        text=True,
        capture_output=True,
        check=True,
    )
    match = re.search(r"^\s*java\.home\s*=\s*(.+)$", probe.stderr, flags=re.MULTILINE)
    if not match:
        raise RuntimeError("Could not determine JAVA_HOME from the installed Java runtime")
    return match.group(1).strip()


def main() -> int:
    parser = argparse.ArgumentParser(description="Run P2Rank on a downloaded PIEZO1 structure.")
    parser.add_argument("--pdb-id", required=True)
    parser.add_argument("--version", default="2.5.1")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--chains", help="Optional comma-separated author chain IDs to analyze")
    args = parser.parse_args()

    pdb_id = args.pdb_id.upper()
    rules = read_json(ROOT / "config" / "consensus_rules.json")
    if pdb_id not in rules["active_pdb_ids"]:
        raise SystemExit(f"{pdb_id} is not an active independently analyzed structure")
    tool_dir = ROOT / "tools" / "p2rank" / f"p2rank_{args.version}"
    launcher = tool_dir / ("prank.bat" if os.name == "nt" else "prank")
    if not launcher.exists():
        raise FileNotFoundError(f"Missing {launcher}; run python scripts/bootstrap_p2rank.py")

    input_path = ROOT / "data" / "raw" / "structures" / f"{pdb_id}.pdb"
    if not input_path.exists():
        raise FileNotFoundError(f"Missing {input_path}; run 01_fetch_structures.py")

    run_dir = ROOT / "runs" / "p2rank" / pdb_id
    output_dir = run_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    command = [
        str(launcher),
        "predict",
        "-f",
        str(input_path),
        "-c",
        "alphafold",
        "-threads",
        str(args.threads),
        "-visualizations",
        "0",
        "-o",
        str(output_dir),
    ]
    if args.chains:
        command.extend(["-chains", args.chains])
    command_text = subprocess.list2cmdline(command)
    metadata = {
        "pdb_id": pdb_id,
        "p2rank_version": args.version,
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": command_text,
        "input_path": str(input_path.relative_to(ROOT)),
        "input_sha256": sha256_file(input_path),
        "reason_for_profile": "P2Rank recommends the alphafold profile for cryo-EM structures because B-factor semantics differ.",
    }
    run_dir.mkdir(parents=True, exist_ok=True)
    write_json(run_dir / "run.json", metadata)

    environment = os.environ.copy()
    environment["JAVA_HOME"] = discover_java_home()
    launch_command = ["cmd.exe", "/d", "/c", command_text] if os.name == "nt" else command
    process = subprocess.run(
        launch_command,
        cwd=tool_dir,
        text=True,
        capture_output=True,
        env=environment,
        check=False,
    )
    (run_dir / "stdout.log").write_text(process.stdout, encoding="utf-8")
    (run_dir / "stderr.log").write_text(process.stderr, encoding="utf-8")
    metadata["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    metadata["returncode"] = process.returncode
    prediction_path = output_dir / f"{pdb_id}.pdb_predictions.csv"
    if prediction_path.exists():
        metadata["prediction_path"] = str(prediction_path.relative_to(ROOT))
        metadata["prediction_sha256"] = sha256_file(prediction_path)
    write_json(run_dir / "run.json", metadata)
    print(process.stdout)
    if process.returncode:
        print(process.stderr)
        raise SystemExit(process.returncode)
    print(f"P2Rank output: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
