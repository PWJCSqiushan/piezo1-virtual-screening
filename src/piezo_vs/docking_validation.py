"""GNINA input/pose reconciliation. Missing or corrupt records cannot be success."""
from __future__ import annotations
import json
import math
from pathlib import Path
from .candidate_selection import read_csv_rows, sha256_file

SCORES = ("minimizedAffinity", "CNNscore", "CNNaffinity")


def validate_prepared_inputs(directory: Path):
    config = json.loads((directory / "docking_config.json").read_text())
    meta = json.loads((directory / "input_run.json").read_text())
    files = {"receptor_pdb": "receptor.pdb", "ligands_sdf": "ligands.sdf",
             "candidate_manifest_csv": "candidate_manifest.csv", "docking_config_json": "docking_config.json"}
    for key, name in files.items():
        if meta.get("output_sha256", {}).get(key) != sha256_file(directory / name):
            raise ValueError("Prepared input hash mismatch: " + name)
    rows = read_csv_rows(directory / "candidate_manifest.csv")
    ids = [r["compound_id"] for r in rows]
    if not rows or len(ids) != len(set(ids)) or ids != config["candidate_ids"] or len(rows) != config["requested_candidate_count"]:
        raise ValueError("Candidate manifest count/identity mismatch")
    if config["purpose"] == "formal_screening":
        if config["receptor_policy"]["review_status"] != "approved" or not config["receptor_policy"]["review_note"]:
            raise ValueError("Unreviewed receptor environment")
        for row in rows:
            if row["approved"] != "true" or row["evaluation_complete"] != "true" or row["review_status"] != "approved":
                raise ValueError("Formal candidate is not evaluated and approved")
    for row in rows:
        if row["pdb_id"] != config["pdb_id"] or row["consensus_id"] != config["consensus_id"]:
            raise ValueError("Candidate structure/site mismatch")
    if any(r["preparation_status"] != "prepared" for r in rows):
        raise ValueError("Partial ligand preparation; inspect candidate_manifest.csv")
    from rdkit import Chem
    molecules = list(Chem.SDMolSupplier(str(directory / "ligands.sdf"), removeHs=False))
    if len(molecules) != len(rows) or any(m is None for m in molecules):
        raise ValueError("Prepared SDF is damaged or incomplete")
    actual_ids = [m.GetProp("_Name") for m in molecules]
    if actual_ids != ids:
        raise ValueError("Prepared SDF identity/order mismatch")
    return config, rows


def inspect_poses(path: Path, candidates: list[dict], run_id: str):
    from rdkit import Chem
    expected = {r["compound_id"]: r for r in candidates}
    valid, issues, seen = [], [], {}
    if path.is_file() and path.stat().st_size:
        supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
    else:
        supplier = []
        issues.append("missing_or_empty_output_sdf")
    visited = 0
    for index, mol in enumerate(supplier, 1):
        visited = index
        try:
            if mol is None or mol.GetNumAtoms() == 0:
                raise ValueError("unreadable_pose")
            cid = mol.GetProp("compound_id") if mol.HasProp("compound_id") else mol.GetProp("_Name")
            if cid not in expected:
                raise ValueError("unexpected_compound_id:" + cid)
            if mol.HasProp("_Name") and mol.GetProp("_Name") != cid:
                raise ValueError("pose_name_identity_conflict:" + cid)
            if not mol.GetNumConformers() or not mol.GetConformer().Is3D():
                raise ValueError("missing_3d_conformer:" + cid)
            xyz = mol.GetConformer().GetPositions()
            if not all(math.isfinite(float(v)) for point in xyz for v in point):
                raise ValueError("nonfinite_pose_coordinate:" + cid)
            Chem.SanitizeMol(mol)
            graph = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
            reference = Chem.MolToSmiles(Chem.MolFromSmiles(expected[cid]["canonical_smiles"]), isomericSmiles=True)
            if graph != reference:
                raise ValueError("pose_structure_changed:" + cid)
            for field in ("pdb_id", "consensus_id", "site_direction"):
                if mol.HasProp(field) and mol.GetProp(field) != expected[cid][field]:
                    raise ValueError("pose_site_identity_conflict:" + cid)
            scores = {key: float(mol.GetProp(key)) for key in SCORES}
            if not all(math.isfinite(v) for v in scores.values()):
                raise ValueError("nonfinite_pose_score:" + cid)
            seen[cid] = seen.get(cid, 0) + 1
            valid.append({**expected[cid], "run_id": run_id, "pose_index": seen[cid], **scores,
                          "interpretation": "computational_pose_not_experimental_binding"})
        except (ValueError, RuntimeError, KeyError) as exc:
            issues.append(f"record_{index}:{exc}")
    if path.is_file():
        raw = path.read_text(encoding="utf-8", errors="replace")
        raw_count = sum(bool(piece.strip()) for piece in raw.split("$$$$"))
        if raw_count != visited or (raw.strip() and not raw.rstrip().endswith("$$$$")):
            issues.append("sdf_record_count_mismatch_or_truncated_tail")
    ledger = [{**r, "run_id": run_id, "valid_pose_count": seen.get(r["compound_id"], 0),
               "docking_status": "computed" if seen.get(r["compound_id"]) else "failed_missing_valid_pose"}
              for r in candidates]
    best = []
    for cid in sorted(seen):
        poses = [r for r in valid if r["compound_id"] == cid]
        best.append(min(poses, key=lambda r: (-r["CNNscore"], r["minimizedAffinity"], r["pose_index"])))
    if any(not r["valid_pose_count"] for r in ledger):
        issues.append("missing_candidate_poses")
    return valid, best, ledger, issues


def check_technical_evidence(paths, config, receptor_hash, binary_hash, batch_ids):
    """Require three successful reviewed molecules and coverage of this environment."""
    ids, environment_covered = set(), False
    for path in paths:
        data = json.loads(path.read_text())
        if data.get("status") != "success" or data.get("stage") != "technical_check":
            raise ValueError("Technical evidence did not pass")
        if data.get("purpose") != config["purpose"] or data["gnina_sha256"] != binary_hash:
            raise ValueError("Technical evidence purpose or binary mismatch")
        for filename, key in (("docked.sdf", "docked_sdf_sha256"), ("candidate_results.csv", "candidate_results_sha256")):
            if sha256_file(path.parent / filename) != data["outputs"][key]:
                raise ValueError("Technical output changed")
        ledger = read_csv_rows(path.parent / "candidate_results.csv")
        observed = {r["compound_id"] for r in ledger if r["docking_status"] == "computed"}
        if not 1 <= len(ledger) <= 3 or observed != set(data["candidate_ids"]) or len(observed) != len(ledger):
            raise ValueError("Technical result ledger is incomplete")
        ids.update(data["candidate_ids"])
        if data["input_sha256"]["receptor_pdb"] == receptor_hash and data["box"] == config["box"] and data["pdb_id"] == config["pdb_id"] and data["consensus_id"] == config["consensus_id"]:
            if not set(data["candidate_ids"]) <= batch_ids:
                raise ValueError("Technical candidates are absent from this selected batch")
            environment_covered = True
    if len(ids) < 3 or not environment_covered:
        raise ValueError("Need at least 3 successful technical-check molecules and matching receptor/site coverage")
