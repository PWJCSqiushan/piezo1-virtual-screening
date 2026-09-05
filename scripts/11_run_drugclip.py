from __future__ import annotations

import argparse
import importlib
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.drugclip_io import inspect_lmdb
from piezo_vs.io_utils import read_json, sha256_file, write_json


def dependency_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    for module_name in ("torch", "numpy", "lmdb", "rdkit", "sklearn", "unicore"):
        try:
            module = importlib.import_module(module_name)
        except ImportError as exc:
            raise RuntimeError(f"Missing DrugCLIP dependency: {module_name}") from exc
        versions[module_name] = str(
            getattr(module, "__version__", getattr(module, "version", "unknown"))
        )
    return versions


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one-pocket DrugCLIP retrieval with provenance.")
    parser.add_argument("--drugclip-root", type=Path, default=ROOT / "third_party" / "DrugCLIP-main")
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--mol-lmdb", type=Path, required=True)
    parser.add_argument("--pocket-lmdb", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--pdb-id", required=True)
    parser.add_argument("--consensus-id", required=True)
    parser.add_argument("--tier", choices=("T1", "T2"), required=True)
    parser.add_argument("--gpu-id", type=int, default=0)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--no-fp16", action="store_true")
    args = parser.parse_args()

    drugclip_root = args.drugclip_root.resolve()
    checkpoint = args.checkpoint.resolve()
    mol_lmdb = args.mol_lmdb.resolve()
    pocket_lmdb = args.pocket_lmdb.resolve()
    output_dir = args.output_dir.resolve()
    retrieval_script = drugclip_root / "unimol" / "retrieval.py"
    dictionary_dir = drugclip_root / "data"
    for required in (
        retrieval_script,
        dictionary_dir / "dict_mol.txt",
        dictionary_dir / "dict_pkt.txt",
        checkpoint,
        mol_lmdb,
        pocket_lmdb,
    ):
        if not required.is_file():
            raise FileNotFoundError(required)
    if checkpoint.stat().st_size < 100_000_000:
        raise ValueError(
            f"Checkpoint is only {checkpoint.stat().st_size} bytes; this is probably a partial download"
        )
    tool_sources = read_json(ROOT / "config" / "tool_sources.json")
    expected_checkpoint_hash = tool_sources["tools"]["drugclip"]["checkpoint_sha256"]
    actual_checkpoint_hash = sha256_file(checkpoint)
    if actual_checkpoint_hash != expected_checkpoint_hash:
        raise ValueError(
            "DrugCLIP checkpoint SHA256 mismatch; refusing to load an unpinned pickle-based checkpoint"
        )
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    embedding_dir = output_dir / "embeddings"
    embedding_dir.mkdir()

    versions = dependency_versions()
    import torch

    if not torch.cuda.is_available():
        raise RuntimeError(
            "This upstream DrugCLIP retrieval path calls move_to_cuda unconditionally; a CUDA GPU is required"
        )
    mol_count, mol_example = inspect_lmdb(mol_lmdb)
    pocket_count, pocket_example = inspect_lmdb(pocket_lmdb)
    if mol_count < 1 or not {"atoms", "coordinates", "smi"} <= set(mol_example):
        raise ValueError("Molecule LMDB is empty or missing atoms/coordinates/smi")
    if pocket_count != 1:
        raise ValueError(
            f"Exactly one pocket is required per retrieval to preserve attribution; found {pocket_count}"
        )
    if not {"pocket", "pocket_atoms", "pocket_coordinates"} <= set(pocket_example):
        raise ValueError("Pocket LMDB is missing pocket/pocket_atoms/pocket_coordinates")
    expected_pocket_name = f"{args.pdb_id.upper()}:{args.consensus_id.upper()}"
    pocket_name = str(pocket_example["pocket"])
    if pocket_name != expected_pocket_name:
        raise ValueError(f"Pocket name {pocket_name!r} does not match {expected_pocket_name!r}")

    command = [
        sys.executable,
        str(retrieval_script),
        "--user-dir",
        str(drugclip_root / "unimol"),
        str(dictionary_dir),
        "--valid-subset",
        "test",
        "--results-path",
        str(output_dir),
        "--num-workers",
        "0",
        "--ddp-backend",
        "c10d",
        "--batch-size",
        str(args.batch_size),
        "--task",
        "drugclip",
        "--loss",
        "in_batch_softmax",
        "--arch",
        "drugclip",
        "--max-pocket-atoms",
        "256",
        "--fp16-init-scale",
        "4",
        "--fp16-scale-window",
        "256",
        "--seed",
        "1",
        "--path",
        str(checkpoint),
        "--log-interval",
        "100",
        "--log-format",
        "simple",
        "--mol-path",
        str(mol_lmdb),
        "--pocket-path",
        str(pocket_lmdb),
        "--emb-dir",
        str(embedding_dir),
    ]
    if not args.no_fp16:
        command.append("--fp16")
    environment = os.environ.copy()
    environment["CUDA_VISIBLE_DEVICES"] = str(args.gpu_id)
    compatibility_dir = ROOT / "compat" / "drugclip"
    previous_pythonpath = environment.get("PYTHONPATH", "")
    environment["PYTHONPATH"] = os.pathsep.join(
        part for part in (str(compatibility_dir), previous_pythonpath) if part
    )
    started_at = datetime.now(timezone.utc)
    metadata = {
        "schema_version": 1,
        "status": "running",
        "started_at_utc": started_at.isoformat(),
        "pdb_id": args.pdb_id.upper(),
        "consensus_id": args.consensus_id.upper(),
        "tier": args.tier,
        "molecule_count": mol_count,
        "pocket_name": pocket_name,
        "pocket_atom_count": len(pocket_example["pocket_atoms"]),
        "versions": versions,
        "gpu_name": torch.cuda.get_device_name(0),
        "checkpoint_loading": {
            "weights_only": True,
            "compatibility_allowlist": str(
                (compatibility_dir / "sitecustomize.py").relative_to(ROOT)
            ),
        },
        "command": command,
        "input_sha256": {
            "checkpoint": actual_checkpoint_hash,
            "mols_lmdb": sha256_file(mol_lmdb),
            "pocket_lmdb": sha256_file(pocket_lmdb),
        },
    }
    metadata_path = output_dir / "run.json"
    write_json(metadata_path, metadata)
    process = subprocess.run(
        command,
        cwd=drugclip_root,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    (output_dir / "stdout.log").write_text(process.stdout, encoding="utf-8")
    (output_dir / "stderr.log").write_text(process.stderr, encoding="utf-8")
    ranked_path = embedding_dir / "ranked_compounds.txt"
    metadata.update(
        {
            "status": "success" if process.returncode == 0 and ranked_path.is_file() else "failed",
            "finished_at_utc": datetime.now(timezone.utc).isoformat(),
            "returncode": process.returncode,
            "ranked_compounds_sha256": sha256_file(ranked_path) if ranked_path.is_file() else "",
        }
    )
    write_json(metadata_path, metadata)
    if metadata["status"] != "success":
        print(process.stderr[-4000:], file=sys.stderr)
        return process.returncode or 1
    print(
        f"DRUGCLIP_RETRIEVAL_OK pdb_id={metadata['pdb_id']} "
        f"consensus_id={metadata['consensus_id']} molecules={mol_count}"
    )
    print(ranked_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
