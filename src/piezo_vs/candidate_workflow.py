"""File-bound Step 8 imports and selection snapshots for finite candidate pools."""
from __future__ import annotations

import json
import math
import shutil
from collections import Counter
from pathlib import Path

from .candidate_selection import (
    TOOL_NAMES, REVIEW_TEMPLATE_FIELDS, CandidateValidationError,
    assess_candidate_rows, read_csv_rows, review_template_rows,
    select_candidates, sha256_file, write_csv_rows, write_external_input_bundle,
)


def import_predictions(manifest_path: Path, master_path: Path, archive_dir: Path):
    """Import declared raw CSV columns, retaining all raw bytes and columns.

    A manifest is required because web exporters do not share a stable schema.
    The operator declares the actual ID/SMILES and prediction column names.
    Files with missing predictions remain failed; a status label alone cannot pass.
    """
    manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    tool = manifest["tool"]
    if tool not in TOOL_NAMES:
        raise CandidateValidationError("Unknown prediction tool")
    if manifest.get("input_master_sha256") != sha256_file(master_path):
        raise CandidateValidationError("External results belong to a different candidate master")
    for key in ("version", "retrieved_at_utc", "files"):
        if not manifest.get(key):
            raise CandidateValidationError("External manifest missing " + key)
    master = {r["compound_id"]: r for r in read_csv_rows(master_path)}
    imported, seen = [], set()
    for spec in manifest["files"]:
        path = (manifest_path.parent / spec["path"]).resolve()
        if not path.is_file() or sha256_file(path) != spec["sha256"]:
            raise CandidateValidationError("Raw prediction file hash mismatch: " + str(path))
        fields = spec["prediction_columns"]
        if not fields or set(fields) & {spec["id_column"], spec.get("smiles_column"), spec.get("status_column")}:
            raise CandidateValidationError("Declare actual prediction columns, distinct from identity/status")
        types = spec.get("prediction_types", {})
        if set(types) != set(fields) or any(t not in {"number", "category"} for t in types.values()):
            raise CandidateValidationError("Declare number/category types for every prediction column")
        rows = read_csv_rows(path)
        status_fields = [k for k in (rows[0] if rows else {}) if k.strip().lower() in {"status", "run_status", "prediction_status", "result_status"}]
        if status_fields and not spec.get("status_column"):
            raise CandidateValidationError("Raw export contains status columns; declare status_column and success_values")
        if not rows or any(k not in rows[0] for k in [spec["id_column"], *fields]):
            raise CandidateValidationError("Raw export lacks declared columns or rows")
        for index, raw in enumerate(rows, 1):
            cid = raw[spec["id_column"]].strip()
            if cid not in master or cid in seen:
                raise CandidateValidationError("Unknown or duplicate external compound: " + cid)
            seen.add(cid)
            if spec.get("smiles_column"):
                from rdkit import Chem
                mol = Chem.MolFromSmiles(raw[spec["smiles_column"]])
                if mol is None or Chem.MolToSmiles(mol, isomericSmiles=True) != master[cid]["canonical_smiles"]:
                    raise CandidateValidationError("External structure change requires reconciliation: " + cid)
            payload = {key: raw[key] for key in fields}
            valid = all(str(v).strip() for v in payload.values())
            for key, value in payload.items():
                if str(value).strip().lower() in {"n/a", "na", "none", "null", "error", "failed", "pending", "-", "--"}:
                    valid = False
                if types[key] == "number":
                    try:
                        number = float(value)
                    except (ValueError, TypeError):
                        valid = False
                    else:
                        if not math.isfinite(number):
                            raise CandidateValidationError("Nonfinite external prediction: " + cid)
                else:
                    allowed = spec.get("allowed_values", {}).get(key, [])
                    if not allowed:
                        raise CandidateValidationError("Category prediction column needs allowed_values: " + key)
                    valid = valid and value in allowed
            if spec.get("status_column"):
                for key in set(status_fields + [spec["status_column"]]):
                    value = raw.get(key, "")
                    valid = valid and value in spec["success_values"] and value.strip().lower() not in {"pending","failed","error","not_started","running","queued"}
            imported.append({
                "compound_id": cid, "canonical_smiles": master[cid]["canonical_smiles"],
                "result_status": "completed" if valid else "failed",
                "prediction_payload_json": json.dumps(payload, ensure_ascii=False, sort_keys=True),
                "all_raw_columns_json": json.dumps(raw, ensure_ascii=False, sort_keys=True),
                "raw_file_sha256": spec["sha256"], "raw_row_number": index,
                "version": manifest["version"], "retrieved_at_utc": manifest["retrieved_at_utc"],
                "structure_check": "matched" if spec.get("smiles_column") else "id_only_no_returned_structure",
            })
        archive_dir.mkdir(parents=True, exist_ok=True)
        target = archive_dir / (spec["sha256"] + path.suffix)
        if not target.exists():
            shutil.copyfile(path, target)
    shutil.copyfile(manifest_path, archive_dir / (sha256_file(manifest_path) + ".manifest.json"))
    return tool, imported


def verified_ranking(path: Path, master_path: Path, pdb_id: str, consensus_id: str):
    metadata = json.loads(path.with_suffix(path.suffix + ".run.json").read_text())
    for key, expected in {"pdb_id": pdb_id, "consensus_id": consensus_id, "tier": "T2"}.items():
        if metadata.get(key) != expected:
            raise CandidateValidationError("Local ranking metadata site/tier mismatch: " + key)
    if metadata["output_sha256"] != sha256_file(path):
        raise CandidateValidationError("Local ranking hash mismatch")
    input_path = Path(metadata["input_run_path"])
    if sha256_file(input_path) != metadata["input_run_sha256"]:
        raise CandidateValidationError("Local ranking preparation changed")
    inputs = json.loads(input_path.read_text())
    if inputs["source_files"]["compounds_sha256"] != sha256_file(master_path):
        raise CandidateValidationError("Local ranking candidate pool mismatch")
    run = json.loads((path.parent / "run.json").read_text())
    rows = read_csv_rows(path)
    if run.get("status") != "success" or run.get("molecule_count") != len(rows) or metadata["row_count"] != len(rows):
        raise CandidateValidationError("Ranking runtime failed or count mismatch")
    for data in (inputs, run, *rows):
        if any(data.get(k) != v for k,v in {"pdb_id":pdb_id,"consensus_id":consensus_id,"tier":"T2"}.items()):
            raise CandidateValidationError("Ranking runtime/row site mismatch")
    for filename, key in (("mols.lmdb","mols_lmdb"),("pocket.lmdb","pocket_lmdb")):
        actual = sha256_file(input_path.parent/filename)
        if actual != inputs["outputs"][key+"_sha256"] or actual != run["input_sha256"][key]:
            raise CandidateValidationError("Ranking LMDB provenance mismatch")
    root = Path(__file__).resolve().parents[2]
    expected_model = json.loads((root/"config/tool_sources.json").read_text())["tools"]["drugclip"]["checkpoint_sha256"]
    if run["input_sha256"]["checkpoint"] != expected_model:
        raise CandidateValidationError("Local ranking used an unpinned model")
    for label in ("structure", "consensus"):
        if sha256_file(root/inputs["source_files"][label]) != inputs["source_files"][label+"_sha256"]:
            raise CandidateValidationError("Ranking receptor/consensus source changed")
    return rows


SELECTION_FIELDS = [
    "compound_id", "canonical_smiles", "site_direction", "pdb_id", "consensus_id",
    "review_status", "approved", "evaluation_complete", "source_rank", "source_drugclip_score",
    "selection_reason", "selection_parameters", "reviewer", "reviewed_at_utc",
    "assessment_comment", "assessment_evidence_sha256", "source_ids", "source_tasks",
]


def build_snapshot(batch: Path, output: Path, external_manifests=(), human_review=None):
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Refusing to overwrite an assessment snapshot")
    output.mkdir(parents=True, exist_ok=True)
    master_path = batch / "audit/candidate_master.csv"
    audit_dir = batch / ("audit_recheck" if (batch / "audit_recheck/audit_summary.json").exists() else "audit")
    audit = json.loads((audit_dir / "audit_summary.json").read_text())
    if sha256_file(batch / "handoff_original.zip") != audit["archive_sha256"] or sha256_file(master_path) != sha256_file(audit_dir / "candidate_master_audited.csv"):
        raise CandidateValidationError("Archived handoff or audited structure inventory changed")
    if audit["exploratory_pool_status"] != "eligible" or sha256_file(master_path) != sha256_file(batch / "audit/candidate_master_audited.csv"):
        raise CandidateValidationError("Candidate audit did not pass or master changed")
    rankings = {
        "C006_C007": verified_ranking(batch / "rerank/8YEZ_C006_two_plus/ranked_compounds.csv", master_path, "8YEZ", "C006"),
        "C016": verified_ranking(batch / "rerank/8ZU3_C016_two_plus/ranked_compounds.csv", master_path, "8ZU3", "C016"),
    }
    external = {}
    for path in external_manifests:
        tool, rows = import_predictions(path, master_path, output / "raw_imports")
        external.setdefault(tool, []).extend(rows)
    master = read_csv_rows(master_path)
    # A ranking ID must still designate the same exact stereochemical structure.
    identities = {r["compound_id"]: r["canonical_smiles"] for r in master}
    for rows in rankings.values():
        for row in rows:
            if identities.get(row["compound_id"]) != row["canonical_smiles"]:
                raise CandidateValidationError("Local ranking identity conflict")
    assessed = assess_candidate_rows(master, rankings, external_rows_by_tool=external,
                                    review_rows=read_csv_rows(human_review) if human_review else None)
    write_csv_rows(output / "candidate_assessment.csv", assessed)
    write_csv_rows(output / "human_review_template.csv", review_template_rows(assessed), REVIEW_TEMPLATE_FIELDS)
    write_external_input_bundle(master, output / "external_inputs")
    selected, selection = select_candidates(assessed)
    for row in selected:
        direction = row["selection_direction"]
        row.update(site_direction=direction, pdb_id="8YEZ" if direction == "C006_C007" else "8ZU3",
                   consensus_id="C006" if direction == "C006_C007" else "C016",
                   review_status="approved", approved=True, evaluation_complete=True)
    write_csv_rows(output / "selected_candidates.csv", selected, SELECTION_FIELDS)
    for pdb, site in (("8YEZ", "C006"), ("8ZU3", "C016")):
        write_csv_rows(output / f"selected_{pdb}_{site}.csv", [r for r in selected if r["pdb_id"] == pdb], SELECTION_FIELDS)
    # Three reviewed molecules cover both intended receptor/site environments.
    technical=[]
    for pdb, limit in (("8YEZ",2),("8ZU3",1)):
        technical.extend(sorted([r for r in selected if r["pdb_id"]==pdb],key=lambda r:r["source_rank"])[:limit])
    for pdb, site in (("8YEZ","C006"),("8ZU3","C016")):
        write_csv_rows(output / f"technical_{pdb}_{site}.csv",[r for r in technical if r["pdb_id"]==pdb],SELECTION_FIELDS)
    # Deterministic replay is checked on the actual snapshot, including a partial/empty result.
    replay, replay_summary = select_candidates(list(reversed(assessed)))
    if [r["compound_id"] for r in selected] != [r["compound_id"] for r in replay] or selection != replay_summary:
        raise AssertionError("Selection replay changed with source row order")
    selection["selected_files"] = {p.name: sha256_file(p) for p in output.glob("*.csv") if p.name.startswith(("selected_", "technical_"))}
    selection.update(selection_sha256=sha256_file(output / "selected_candidates.csv"),
                     assessment_sha256=sha256_file(output / "candidate_assessment.csv"),
                     deterministic_replay_passed=True)
    (output / "selection_manifest.json").write_text(json.dumps(selection, indent=2), encoding="utf-8")
    status = {
        "candidate_count": len(master), "local_descriptors": dict(Counter(r["rdkit_status"] for r in assessed)),
        "external_tools": {t: dict(Counter(r[t+"_result_state"] for r in assessed)) for t in TOOL_NAMES},
        "human_review": dict(Counter(r["human_approval_status"] for r in assessed)),
        "selected_count": len(selected), "requested_count": 20, "selection_status": selection["status"],
        "gnina_technical_check": "not_started", "gnina_formal": "not_started",
        "failed_local_descriptors": sum(r["rdkit_status"] != "parsed" for r in assessed),
        "source_master_sha256": sha256_file(master_path),
        "method": "finite returned pool; local reranking; predictions and human review precede docking",
    }
    (output / "status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    return status
