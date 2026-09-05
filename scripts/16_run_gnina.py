from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.docking import numeric_property
from piezo_vs.io_utils import read_json, sha256_file, write_json


POSE_FIELDS = [
    "compound_id",
    "source",
    "source_id",
    "drugclip_rank",
    "drugclip_score",
    "pose_index",
    "minimizedAffinity",
    "CNNscore",
    "CNNaffinity",
    "CNN_VS",
    "CNNaffinity_variance",
    "interpretation",
]


def molecule_properties(molecule) -> dict[str, str]:
    properties = {name: molecule.GetProp(name) for name in molecule.GetPropNames()}
    properties.setdefault("compound_id", molecule.GetProp("_Name") if molecule.HasProp("_Name") else "")
    return properties


def parse_docked_sdf(path: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise RuntimeError("GNINA result parsing requires RDKit") from exc
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
    pose_counts: dict[str, int] = {}
    rows: list[dict[str, object]] = []
    for index, molecule in enumerate(supplier, 1):
        if molecule is None:
            continue
        properties = molecule_properties(molecule)
        compound_id = properties.get("compound_id", "")
        if not compound_id:
            compound_id = f"unknown_pose_{index}"
        pose_counts[compound_id] = pose_counts.get(compound_id, 0) + 1
        row: dict[str, object] = {
            "compound_id": compound_id,
            "source": properties.get("source", ""),
            "source_id": properties.get("source_id", ""),
            "drugclip_rank": properties.get("rank", properties.get("drugclip_rank", "")),
            "drugclip_score": properties.get("drugclip_score", ""),
            "pose_index": pose_counts[compound_id],
            "interpretation": "computational_docking_pose_not_experimental_binding_evidence",
        }
        for field in (
            "minimizedAffinity",
            "CNNscore",
            "CNNaffinity",
            "CNN_VS",
            "CNNaffinity_variance",
        ):
            row[field] = properties.get(field, "")
        rows.append(row)
    if not rows:
        raise ValueError("GNINA output SDF contains no readable poses")

    def sort_key(row: dict[str, object]) -> tuple[float, float, int]:
        cnn_score = numeric_property(row, "CNNscore")
        affinity = numeric_property(row, "minimizedAffinity")
        return (
            -(cnn_score if cnn_score is not None else float("-inf")),
            affinity if affinity is not None else float("inf"),
            int(row["pose_index"]),
        )

    best = []
    for compound_id in sorted(pose_counts):
        candidates = [row for row in rows if row["compound_id"] == compound_id]
        best.append(sorted(candidates, key=sort_key)[0])
    best.sort(key=sort_key)
    for rank, row in enumerate(best, 1):
        row["docking_rank"] = rank
    return rows, best


def write_rows(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run pinned GNINA and export traceable pose tables.")
    parser.add_argument("--gnina", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--num-modes", type=int, default=9)
    parser.add_argument("--cpu", type=int, default=8)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--gpu-id", type=int, default=0)
    parser.add_argument("--no-gpu", action="store_true")
    args = parser.parse_args()
    if args.exhaustiveness < 1 or args.num_modes < 1 or args.cpu < 1:
        raise SystemExit("exhaustiveness, num-modes and cpu must be positive")
    gnina = args.gnina.resolve()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    receptor = input_dir / "receptor.pdb"
    ligands = input_dir / "ligands.sdf"
    config_path = input_dir / "docking_config.json"
    input_run_path = input_dir / "input_run.json"
    for path in (gnina, receptor, ligands, config_path, input_run_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    tool_config = read_json(ROOT / "config" / "tool_sources.json")["tools"]["gnina"]
    actual_hash = sha256_file(gnina)
    if gnina.stat().st_size != int(tool_config["binary_size_bytes"]):
        raise ValueError("GNINA binary size mismatch; refusing to execute a partial download")
    if actual_hash != tool_config["binary_sha256"]:
        raise ValueError("GNINA binary SHA256 mismatch; refusing to execute an unpinned binary")
    config = read_json(config_path)
    box = config["box"]
    docked_sdf = output_dir / "docked.sdf"
    command = [
        str(gnina),
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligands),
        "--center_x",
        str(box["center_x"]),
        "--center_y",
        str(box["center_y"]),
        "--center_z",
        str(box["center_z"]),
        "--size_x",
        str(box["size_x"]),
        "--size_y",
        str(box["size_y"]),
        "--size_z",
        str(box["size_z"]),
        "--out",
        str(docked_sdf),
        "--cnn_scoring",
        "rescore",
        "--pose_sort_order",
        "CNNscore",
        "--exhaustiveness",
        str(args.exhaustiveness),
        "--num_modes",
        str(args.num_modes),
        "--cpu",
        str(args.cpu),
        "--seed",
        str(args.seed),
    ]
    if args.no_gpu:
        command.append("--no_gpu")
    else:
        command.extend(["--device", str(args.gpu_id)])
    version_process = subprocess.run([str(gnina), "--version"], text=True, capture_output=True, check=False)
    started_at = datetime.now(timezone.utc)
    run_metadata = {
        "schema_version": 1,
        "status": "running",
        "started_at_utc": started_at.isoformat(),
        "pdb_id": config["pdb_id"],
        "consensus_id": config["consensus_id"],
        "tier": config["tier"],
        "support_count": config["support_count"],
        "gnina_version": (version_process.stdout or version_process.stderr).strip(),
        "gnina_sha256": actual_hash,
        "command": command,
        "input_sha256": {
            "receptor_pdb": sha256_file(receptor),
            "ligands_sdf": sha256_file(ligands),
            "docking_config": sha256_file(config_path),
            "input_run": sha256_file(input_run_path),
        },
        "scientific_status": "computational_docking_prediction_only",
    }
    metadata_path = output_dir / "run.json"
    write_json(metadata_path, run_metadata)
    process = subprocess.run(command, text=True, capture_output=True, check=False)
    (output_dir / "stdout.log").write_text(process.stdout, encoding="utf-8")
    (output_dir / "stderr.log").write_text(process.stderr, encoding="utf-8")
    run_metadata.update(
        {
            "finished_at_utc": datetime.now(timezone.utc).isoformat(),
            "returncode": process.returncode,
            "status": "success" if process.returncode == 0 and docked_sdf.is_file() else "failed",
        }
    )
    if run_metadata["status"] != "success":
        write_json(metadata_path, run_metadata)
        print(process.stderr[-4000:], file=sys.stderr)
        return process.returncode or 1

    poses, best = parse_docked_sdf(docked_sdf)
    all_poses_path = output_dir / "all_poses.csv"
    best_path = output_dir / "best_poses.csv"
    write_rows(all_poses_path, poses, POSE_FIELDS)
    write_rows(best_path, best, ["docking_rank", *POSE_FIELDS])
    run_metadata["outputs"] = {
        "docked_sdf_sha256": sha256_file(docked_sdf),
        "all_poses_csv_sha256": sha256_file(all_poses_path),
        "best_poses_csv_sha256": sha256_file(best_path),
        "pose_count": len(poses),
        "compound_count": len(best),
        "best_pose_rule": "highest CNNscore then lowest minimizedAffinity within each compound",
    }
    write_json(metadata_path, run_metadata)
    print(
        f"GNINA_DOCKING_OK pdb_id={config['pdb_id']} consensus_id={config['consensus_id']} "
        f"compounds={len(best)} poses={len(poses)}"
    )
    print(best_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
