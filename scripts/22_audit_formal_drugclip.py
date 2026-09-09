from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors


ALLOWED_ORGANIC_ELEMENTS = {"H", "B", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}
TASK_FILES = (
    "input_run.json",
    "molecule_manifest.csv",
    "ranked_compounds.csv",
    "run.log",
    "standardized_pocket.pdb",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def pocket_atom_hash(path: Path) -> tuple[str, int]:
    atom_lines = []
    with path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")):
                atom_lines.append(line.rstrip("\r\n"))
    payload = ("\n".join(atom_lines) + "\n").encode("ascii", errors="replace")
    return hashlib.sha256(payload).hexdigest().upper(), len(atom_lines)


def molecule_flags(smiles: str) -> tuple[list[str], list[str], dict[str, object]]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ["rdkit_parse_failed"], [], {}
    atoms = list(mol.GetAtoms())
    elements = {atom.GetSymbol() for atom in atoms}
    heavy_atoms = mol.GetNumHeavyAtoms()
    carbon_atoms = sum(atom.GetSymbol() == "C" for atom in atoms)
    fragments = len(Chem.GetMolFrags(mol))
    isotope_atoms = sum(atom.GetIsotope() != 0 for atom in atoms)
    review: list[str] = []
    if carbon_atoms == 0:
        review.append("no_carbon")
    if heavy_atoms < 6:
        review.append("very_small_less_than_6_heavy_atoms")
    if fragments > 1:
        review.append("multiple_fragments_or_salt")
    disallowed = sorted(elements - ALLOWED_ORGANIC_ELEMENTS)
    if disallowed:
        review.append("metal_or_nonstandard_element:" + ",".join(disallowed))
    if isotope_atoms:
        review.append("isotope_labeled")
    mw = float(Descriptors.MolWt(mol))
    if mw < 100:
        review.append("mw_below_100")
    if mw > 800:
        review.append("mw_above_800")
    return [], review, {
        "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
        "heavy_atoms": heavy_atoms,
        "carbon_atoms": carbon_atoms,
        "fragments": fragments,
        "elements": ",".join(sorted(elements)),
        "mw": round(mw, 3),
    }


def parse_hash_manifest(path: Path) -> dict[str, str]:
    expected: dict[str, str] = {}
    with path.open(encoding="utf-8-sig") as handle:
        for raw in handle:
            raw = raw.strip()
            if not raw:
                continue
            digest, relative = raw.split(None, 1)
            expected[relative.replace("\\", "/")] = digest.upper()
    return expected


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Audit a teammate DrugCLIP export without treating ranking as biological evidence."
    )
    parser.add_argument("base", type=Path)
    parser.add_argument("repo", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--score-kind",
        choices=("unknown", "similarity", "zscore"),
        default="unknown",
        help="Only zscore enables the optional numeric threshold; unknown is safest.",
    )
    parser.add_argument("--zscore-threshold", type=float, default=3.0)
    args = parser.parse_args()
    base = args.base.resolve()
    repo = args.repo.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir
        else repo / "runs" / "formal_drugclip_audit" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    )
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    for required in ("SHA256SUMS.txt", "candidate_selection.csv", "screening_summary.csv", "README.md"):
        if not (base / required).is_file():
            raise FileNotFoundError(base / required)

    expected_hashes = parse_hash_manifest(base / "SHA256SUMS.txt")
    hash_failures: list[str] = []
    missing_hash_targets: list[str] = []
    for relative, expected in expected_hashes.items():
        path = base / relative
        if not path.is_file():
            missing_hash_targets.append(relative)
        elif sha256(path) != expected:
            hash_failures.append(relative)

    selection = {
        (row["pdb_id"], row["consensus_id"]): row
        for row in read_csv(base / "candidate_selection.csv")
    }
    summary = {
        (row["pdb_id"], row["consensus_id"]): row
        for row in read_csv(base / "screening_summary.csv")
    }

    tasks: dict[str, dict[str, object]] = {}
    task_rows: list[dict[str, object]] = []
    review_rows: list[dict[str, object]] = []
    ranked_sets: dict[str, set[str]] = {}
    ranked_file_hashes: dict[str, list[str]] = {}
    atom_groups: dict[str, list[str]] = {}
    normalized_dir = output_dir / "normalized"

    for task_dir in sorted(path for path in base.iterdir() if path.is_dir()):
        task = task_dir.name
        if "_" not in task:
            continue
        pdb_id, consensus_id = task.split("_", 1)
        missing_files = [name for name in TASK_FILES if not (task_dir / name).is_file()]
        if missing_files:
            tasks[task] = {"hard_blockers": ["missing_task_files:" + ",".join(missing_files)]}
            continue

        ranked_path = task_dir / "ranked_compounds.csv"
        ranked = read_csv(ranked_path)
        manifest = read_csv(task_dir / "molecule_manifest.csv")
        metadata = json.loads((task_dir / "input_run.json").read_text(encoding="utf-8-sig"))
        consensus_path = repo / "results" / "consensus" / f"{pdb_id}_consensus_pockets.csv"
        consensus_rows = read_csv(consensus_path)
        consensus = next((row for row in consensus_rows if row["consensus_id"] == consensus_id), None)
        ids = [row.get("compound_id", "") for row in ranked]
        ranks = [int(row["rank"]) for row in ranked]
        scores = [float(row["drugclip_score"]) for row in ranked]
        ranked_sets[task] = set(ids)

        hard_counts: Counter[str] = Counter()
        review_counts: Counter[str] = Counter()
        canonical_counts: Counter[str] = Counter()
        normalized_rows: list[dict[str, object]] = []
        manifest_by_id = {row.get("compound_id", ""): row for row in manifest}
        for row in ranked:
            compound_id = row.get("compound_id", "")
            hard, review, props = molecule_flags(row.get("smiles", ""))
            canonical = str(props.get("canonical_smiles", ""))
            if canonical:
                canonical_counts[canonical] += 1
            hard_counts.update(hard)
            review_counts.update(review)
            source = manifest_by_id.get(compound_id, {})
            normalized = {
                "pdb_id": pdb_id,
                "consensus_id": consensus_id,
                "tier": metadata.get("tier", ""),
                "rank": int(row["rank"]),
                "compound_id": compound_id,
                "canonical_smiles": canonical or row.get("smiles", ""),
                "source": source.get("source", ""),
                "source_id": source.get("source_id", ""),
                "source_url": source.get("source_url", ""),
                "drugclip_score": float(row["drugclip_score"]),
                "score_kind": args.score_kind,
                "hard_blockers": ";".join(hard),
                "review_flags": ";".join(review),
                "interpretation": "computational_rank_only_not_binding_efficacy_or_safety",
            }
            normalized_rows.append(normalized)
            review_rows.append({**normalized, **props})

        duplicate_canonical_count = sum(
            count - 1 for count in canonical_counts.values() if count > 1
        )
        manifest_pairs = [(r.get("compound_id", ""), r.get("smiles", "")) for r in manifest]
        ranked_pairs = [(r.get("compound_id", ""), r.get("smiles", "")) for r in ranked]
        atom_hash, atom_count = pocket_atom_hash(task_dir / "standardized_pocket.pdb")
        atom_groups.setdefault(atom_hash, []).append(task)
        ranked_hash = sha256(ranked_path)
        ranked_file_hashes.setdefault(ranked_hash, []).append(task)
        task_key = (pdb_id, consensus_id)
        metadata_missing = sorted(
            key for key, value in metadata.get("inputs", {}).items() if not value
        )
        hard_blockers: list[str] = []
        if hard_counts:
            hard_blockers.append("rdkit_parse_failures")
        if metadata_missing:
            hard_blockers.append("missing_input_hashes")
        if task_key not in summary:
            hard_blockers.append("missing_summary_row")
        if task_key not in selection:
            hard_blockers.append("missing_selection_row")
        library = metadata.get("library", {})
        if not library.get("source_url") or not library.get("version"):
            hard_blockers.append("library_provenance_incomplete")
        drugclip = metadata.get("drugclip", {})
        if not drugclip.get("source_commit") or drugclip.get("seed") is None:
            hard_blockers.append("model_run_provenance_incomplete")

        summary_matches_count = (
            task_key in summary and int(summary[task_key]["molecule_count"]) == len(ranked)
        )
        summary_matches_top = (
            task_key in summary
            and bool(ids)
            and summary[task_key]["top_compound"] == ids[0]
            and math.isclose(float(summary[task_key]["top_score"]), scores[0])
        )
        selection_matches = bool(consensus) and task_key in selection and all(
            str(selection[task_key][key]) == str(consensus[key])
            for key in ("tier", "support_count", "pocket_scope")
        )

        conditional_zscore_count = (
            sum(score > args.zscore_threshold for score in scores)
            if args.score_kind == "zscore"
            else None
        )
        info = {
            "row_count": len(ranked),
            "manifest_row_count": len(manifest),
            "ranks_are_exact_1_to_n": ranks == list(range(1, len(ranked) + 1)),
            "scores_nonincreasing": all(a >= b for a, b in zip(scores, scores[1:])),
            "unique_compound_ids": len(ids) == len(set(ids)),
            "manifest_exactly_matches_ranked": manifest_pairs == ranked_pairs,
            "metadata_matches_task": (
                metadata.get("pdb_id") == pdb_id
                and metadata.get("consensus_id") == consensus_id
            ),
            "metadata_input_hashes_missing": metadata_missing,
            "summary_matches_count": summary_matches_count,
            "summary_matches_top": summary_matches_top,
            "selection_matches_authoritative_consensus": selection_matches,
            "pocket_atom_sha256": atom_hash,
            "pocket_atom_count": atom_count,
            "ranked_file_sha256": ranked_hash,
            "score_min": min(scores) if scores else None,
            "score_max": max(scores) if scores else None,
            "score_kind": args.score_kind,
            "zscore_threshold_applied": args.score_kind == "zscore",
            "conditional_zscore_above_threshold_count": conditional_zscore_count,
            "hard_blocker_molecule_count": sum(hard_counts.values()),
            "review_flagged_molecule_count": sum(
                1 for row in normalized_rows if row["review_flags"]
            ),
            "hard_blocker_counts": dict(sorted(hard_counts.items())),
            "review_flag_counts": dict(sorted(review_counts.items())),
            "duplicate_canonical_smiles_count": duplicate_canonical_count,
            "hard_blockers": hard_blockers,
            "current_gate": "blocked_before_gnina" if hard_blockers else "requires_manual_review",
        }
        tasks[task] = info
        task_rows.append({"task": task, **info})
        write_csv(normalized_dir / f"{task}_ranked_normalized.csv", normalized_rows)

    overlaps: dict[str, dict[str, object]] = {}
    task_names = sorted(ranked_sets)
    for index, left in enumerate(task_names):
        for right in task_names[index + 1 :]:
            intersection = len(ranked_sets[left] & ranked_sets[right])
            union = len(ranked_sets[left] | ranked_sets[right])
            overlaps[f"{left}__{right}"] = {
                "intersection": intersection,
                "union": union,
                "jaccard": round(intersection / union, 6) if union else 0.0,
                "comparable": ranked_sets[left] == ranked_sets[right],
                "reason": (
                    "same_compound_universe"
                    if ranked_sets[left] == ranked_sets[right]
                    else "different_compound_universe_do_not_compare_ranks"
                ),
            }

    identical_atom_groups = [members for members in atom_groups.values() if len(members) > 1]
    identical_rank_groups = [members for members in ranked_file_hashes.values() if len(members) > 1]
    result = {
        "schema_version": 2,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "scientific_status": "computational_prediction_only",
        "base": str(base),
        "output_dir": str(output_dir),
        "overall_gate": "blocked_before_gnina",
        "overall_blockers": [
            "compound_library_provenance_not_independently_reproducible",
            "web_score_definition_not_confirmed_as_zscore",
            "8YEZ_C006_and_C007_used_different_compound_universes",
            "C016_inputs_and_rankings_are_identical_across_three_structure_labels",
        ],
        "score_semantics": {
            "declared_kind": args.score_kind,
            "metadata_sort_by": "similarity",
            "zscore_threshold": args.zscore_threshold,
            "threshold_applied": args.score_kind == "zscore",
            "note": (
                "A numeric >3 count is only valid when the exported score is explicitly "
                "confirmed to be a z-score. Current web exports are not so confirmed."
            ),
        },
        "hash_manifest": {
            "entries": len(expected_hashes),
            "missing_targets": missing_hash_targets,
            "mismatches": hash_failures,
            "screening_summary_covered": "screening_summary.csv" in expected_hashes,
        },
        "root_readme_claims_no_ranking": (
            "No DrugCLIP ranking was generated"
            in (base / "README.md").read_text(encoding="utf-8-sig")
        ),
        "tasks": tasks,
        "ranked_set_overlaps": overlaps,
        "identical_pocket_atom_groups": identical_atom_groups,
        "identical_ranked_file_groups": identical_rank_groups,
        "retention_policy": (
            "All rows are retained. Flags request review and are not automatic efficacy, "
            "toxicity, oral-delivery, inhalation-delivery, or exclusion decisions."
        ),
    }
    (output_dir / "audit_summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    task_fields = [
        "task", "row_count", "manifest_row_count", "ranks_are_exact_1_to_n",
        "scores_nonincreasing", "unique_compound_ids", "manifest_exactly_matches_ranked",
        "metadata_matches_task", "summary_matches_count", "summary_matches_top",
        "selection_matches_authoritative_consensus", "score_kind",
        "zscore_threshold_applied", "conditional_zscore_above_threshold_count",
        "hard_blocker_molecule_count", "review_flagged_molecule_count",
        "duplicate_canonical_smiles_count", "current_gate",
    ]
    write_csv(output_dir / "task_audit.csv", task_rows, task_fields)
    review_fields = [
        "pdb_id", "consensus_id", "tier", "rank", "compound_id", "canonical_smiles",
        "source", "source_id", "source_url", "drugclip_score", "score_kind",
        "hard_blockers", "review_flags", "heavy_atoms", "carbon_atoms", "fragments",
        "elements", "mw", "interpretation",
    ]
    write_csv(output_dir / "molecule_review_flags.csv", review_rows, review_fields)

    rerun_rows = [
        {
            "priority": 1,
            "batch": "8YEZ_symmetry_diagnostic",
            "pdb_id": "8YEZ",
            "consensus_id": consensus_id,
            "library_rule": "same_exact_fixed_library_for_both_tasks",
            "purpose": "technical_validation",
            "gate_after_run": "same_universe_rank_and_top_k_comparison",
        }
        for consensus_id in ("C006", "C007")
    ]
    rerun_rows.extend(
        {
            "priority": 2,
            "batch": "C016_identity_resolution",
            "pdb_id": pdb_id,
            "consensus_id": "C016",
            "library_rule": "do_not_rerun_until_structure_identity_and_library_provenance_are_resolved",
            "purpose": "blocked",
            "gate_after_run": "independent_input_hashes_required",
        }
        for pdb_id in ("8ZU3", "8YFC", "9VMX")
    )
    write_csv(output_dir / "rerun_tasks.csv", rerun_rows)

    evidence_request = """# DrugCLIP 网页任务补证清单

请原运行同学按任务逐一补交：

1. 任务页链接或完整截图：任务 ID、上传文件名、所选分子库、参数、提交时间、完成时间。
2. 网页直接下载的原始文件或 API 响应，不要只提供二次整理 CSV。
3. comp 库的正式名称、供应商、版本/日期、全库规模、许可和可获得的原始库文件/哈希。
4. 明确 100/1000 是全库规模、返回条数还是页面导出上限。
5. 明确网页 drugclip_score 的数学定义；在确认是 z-score 前，不使用 >3 作筛选门槛。
6. 说明网页上传的是完整蛋白还是提取口袋，以及是否做了额外标准化。

这些材料用于可复现性门禁，不是否定已完成的网页操作。
"""
    (output_dir / "evidence_request.md").write_text(evidence_request, encoding="utf-8")
    print(
        f"FORMAL_DRUGCLIP_AUDIT_OK tasks={len(tasks)} rows={len(review_rows)} "
        f"gate={result['overall_gate']} score_kind={args.score_kind}"
    )
    print(output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())