from __future__ import annotations

import csv
import html
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from piezo_vs.io_utils import read_json, sha256_file, write_json


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def count_csv_rows(path: Path) -> int:
    return len(read_csv_rows(path)) if path.is_file() else 0


def _safe_read_json(path: Path, default: Any) -> Any:
    try:
        return read_json(path)
    except (FileNotFoundError, json.JSONDecodeError, TypeError):
        return default


def _tool_counts(root: Path, pdb_id: str) -> dict[str, int]:
    p2rank = count_csv_rows(root / "runs" / "p2rank" / pdb_id / "output" / f"{pdb_id}.pdb_predictions.csv")
    dogsite_latest = _safe_read_json(root / "runs" / "dogsite3" / pdb_id / "latest_job.json", {})
    dogsite_result = {}
    if dogsite_latest.get("job_id"):
        dogsite_result = _safe_read_json(
            root
            / "runs"
            / "dogsite3"
            / pdb_id
            / str(dogsite_latest["job_id"])
            / "result.json",
            {},
        )
    dogsite = len(dogsite_result.get("residues", []))
    fpocket_latest = _safe_read_json(root / "runs" / "fpocket" / pdb_id / "latest_run.json", {})
    fpocket = int(fpocket_latest.get("pocket_count", 0) or 0)
    return {"p2rank": p2rank, "dogsite3": dogsite, "fpocket": fpocket}


def _consensus_summary(root: Path, pdb_id: str, tier_map: dict[int, str]) -> dict[str, Any]:
    path = root / "results" / "consensus" / f"{pdb_id}_consensus_pockets.csv"
    rows = read_csv_rows(path) if path.is_file() else []
    tier_counts = {"T1": 0, "T2": 0}
    scope_counts: dict[str, int] = {}
    errors: list[str] = []
    for row in rows:
        try:
            support = int(row.get("support_count", "0"))
        except ValueError:
            support = 0
        tier = row.get("tier", "")
        if tier in tier_counts:
            tier_counts[tier] += 1
        if tier_map.get(support) != tier:
            errors.append(f"{row.get('consensus_id', '?')}: support_count={support}, tier={tier!r}")
        scope = row.get("pocket_scope", "unknown") or "unknown"
        scope_counts[scope] = scope_counts.get(scope, 0) + 1
    return {
        "path": str(path),
        "count": len(rows),
        "tier_counts": tier_counts,
        "scope_counts": scope_counts,
        "tier_mapping_errors": errors,
        "top_regions_by_support": sorted(
            rows,
            key=lambda row: (
                -int(row.get("support_count", "0") or 0),
                -int(row.get("core_residue_count_2plus", "0") or 0),
                float(row.get("max_member_center_distance", "999") or 999),
            ),
        )[:10],
    }


def build_snapshot(
    root: Path,
    *,
    ranked_csv: Path | None = None,
    input_run: Path | None = None,
    docking_run: Path | None = None,
    assessment_run: Path | None = None,
    stability_json: Path | None = None,
    purpose: str = "technical_validation",
) -> dict[str, Any]:
    structures_config = read_json(root / "config" / "structures.json")
    consensus_rules = read_json(root / "config" / "consensus_rules.json")
    active_ids = list(consensus_rules["active_pdb_ids"])
    tier_map = {
        int(value["support_count"]): tier
        for tier, value in consensus_rules["tiers"].items()
    }
    manifest_path = root / "results" / "structure_manifest.csv"
    manifest_rows = read_csv_rows(manifest_path) if manifest_path.is_file() else []
    manifest_by_id = {row["pdb_id"]: row for row in manifest_rows}

    structures = []
    all_tool_runs = True
    all_consensus = True
    tier_errors: list[str] = []
    coordinate_groups: dict[str, list[str]] = {}
    total_consensus = 0
    for pdb_id in active_ids:
        manifest = manifest_by_id.get(pdb_id, {})
        tools = _tool_counts(root, pdb_id)
        consensus = _consensus_summary(root, pdb_id, tier_map)
        all_tool_runs = all_tool_runs and all(value > 0 for value in tools.values())
        all_consensus = all_consensus and consensus["count"] > 0
        total_consensus += consensus["count"]
        tier_errors.extend(f"{pdb_id}/{error}" for error in consensus["tier_mapping_errors"])
        fingerprint = manifest.get("coordinate_atom_sha256", "")
        if fingerprint:
            coordinate_groups.setdefault(fingerprint, []).append(pdb_id)
        structures.append(
            {
                "pdb_id": pdb_id,
                "category": manifest.get("category", ""),
                "mutation_label": manifest.get("mutation_label", ""),
                "resolution_angstrom": manifest.get("resolution_angstrom", ""),
                "has_mdfic": manifest.get("has_mdfic", ""),
                "modeled_fraction": manifest.get("modeled_fraction", ""),
                "mutation_coordinate_status": manifest.get("mutation_coordinate_status", ""),
                "coordinate_atom_sha256": fingerprint,
                "tool_counts": tools,
                "consensus": consensus,
            }
        )

    alerts: list[dict[str, str]] = []
    for fingerprint, pdb_ids in coordinate_groups.items():
        if len(pdb_ids) > 1:
            alerts.append(
                {
                    "severity": "high",
                    "title": "多个条目使用相同 ATOM HETATM 坐标",
                    "detail": (
                        f"{', '.join(pdb_ids)} 的坐标记录指纹相同（{fingerprint[:12]}...）。"
                        "这些结构必须继续独立留档，但重复口袋结果不能解释为突变特异性差异。"
                    ),
                }
            )
    for structure in structures:
        status = structure["mutation_coordinate_status"]
        if structure["category"] == "mutant" and status not in {
            "substitution_modeled_as_VAL",
            "deletion_site_has_atoms_unexpected",
        }:
            alerts.append(
                {
                    "severity": "high",
                    "title": f"{structure['pdb_id']} 突变坐标证据受限",
                    "detail": f"{structure['mutation_label']}: {status}。不能宣称已观察到突变导致的局部口袋变化。",
                }
            )
    if tier_errors:
        alerts.append(
            {
                "severity": "critical",
                "title": "Tier 与 support_count 映射不一致",
                "detail": "; ".join(tier_errors[:8]),
            }
        )

    input_metadata = _safe_read_json(input_run, {}) if input_run else {}
    ranked_rows = read_csv_rows(ranked_csv) if ranked_csv and ranked_csv.is_file() else []
    docking_metadata = _safe_read_json(docking_run, {}) if docking_run else {}
    assessment_metadata = _safe_read_json(assessment_run, {}) if assessment_run else {}
    stability_metadata = _safe_read_json(stability_json, {}) if stability_json else {}
    docking_best_path = docking_run.parent / "best_poses.csv" if docking_run else None
    docking_best_rows = (
        read_csv_rows(docking_best_path)
        if docking_best_path and docking_best_path.is_file()
        else []
    )
    if ranked_rows and input_metadata:
        for field in ("pdb_id", "consensus_id", "tier"):
            values = {row.get(field, "") for row in ranked_rows}
            if values != {str(input_metadata.get(field, ""))}:
                raise ValueError(
                    f"Ranked CSV {field} values {sorted(values)} do not match input_run.json "
                    f"value {input_metadata.get(field)!r}"
                )
        if len(ranked_rows) != int(input_metadata.get("accepted_compound_count", -1)):
            raise ValueError("Ranked row count does not match accepted_compound_count in input_run.json")
        selected_consensus_path = (
            root / "results" / "consensus" / f"{input_metadata.get('pdb_id', '')}_consensus_pockets.csv"
        )
        expected_consensus_hash = input_metadata.get("source_files", {}).get("consensus_sha256", "")
        if not selected_consensus_path.is_file() or sha256_file(selected_consensus_path) != expected_consensus_hash:
            raise ValueError("DrugCLIP input consensus SHA256 does not match the current consensus table")
        selected_rows = [
            row for row in read_csv_rows(selected_consensus_path)
            if row.get("consensus_id") == input_metadata.get("consensus_id")
        ]
        if len(selected_rows) != 1:
            raise ValueError("Selected consensus ID is missing or duplicated in the current consensus table")
        selected_consensus = selected_rows[0]
    else:
        selected_consensus = {}
    if docking_metadata:
        if docking_metadata.get("status") != "success":
            raise ValueError("Docking run supplied to the report is not successful")
        if ranked_rows and docking_metadata.get("pdb_id") != ranked_rows[0].get("pdb_id"):
            raise ValueError("Docking and DrugCLIP PDB IDs do not match")
        if ranked_rows and docking_metadata.get("consensus_id") != ranked_rows[0].get("consensus_id"):
            raise ValueError("Docking and DrugCLIP consensus IDs do not match")
    technical_only = purpose != "formal_screening"
    stages = [
        {"step": 1, "name": "结构与元数据", "status": "complete" if len(manifest_by_id) == len(active_ids) else "incomplete", "evidence": f"{len(manifest_by_id)}/{len(active_ids)} 个正式结构有清单"},
        {"step": 2, "name": "三工具独立口袋预测", "status": "complete" if all_tool_runs else "incomplete", "evidence": "DoGSite3 fpocket P2Rank 均有本地输出" if all_tool_runs else "至少一个结构缺工具输出"},
        {"step": 3, "name": "同结构共识匹配", "status": "complete" if all_consensus else "incomplete", "evidence": f"共 {total_consensus} 个同结构共识区域"},
        {"step": 4, "name": "support_count 分层", "status": "complete" if all_consensus and not tier_errors else "incomplete", "evidence": "T1=2个工具 T2=3个工具"},
        {"step": 5, "name": "标准化口袋提取", "status": "validated" if input_metadata else "implemented_not_run", "evidence": f"{input_metadata.get('residue_count', 0)} 个残基 {input_metadata.get('pocket_atom_count', 0)} 个原子" if input_metadata else "代码已实现 尚未为本报告指定输入运行"},
        {"step": 6, "name": "共识口袋进入 DrugCLIP", "status": "validated" if ranked_rows else "implemented_not_run", "evidence": f"{len(ranked_rows)} 个分子获得带口袋归属的排序" if ranked_rows else "GPU流程已实现 本报告未附排序"},
        {"step": 7, "name": "化合物库快速排名", "status": "technical_validation_only" if ranked_rows and technical_only else ("complete_run" if ranked_rows else "pending"), "evidence": "工程压力测试 不代表正式候选" if ranked_rows and technical_only else ("指定分子库已完成一次运行" if ranked_rows else "等待正式分子库")},
        {"step": 8, "name": "可开发性与毒性分流", "status": "validated" if assessment_metadata.get("branch_counts") else "implemented_not_run", "evidence": f"已生成 {assessment_metadata.get('compound_count', 0)} 个候选的无损评估入口 等待外部预测" if assessment_metadata else "导入 导出 分流与专家复核队列代码已实现 尚无真实预测导出"},
        {"step": 9, "name": "GNINA docking", "status": "technical_validation_only" if docking_metadata and technical_only else ("complete_run" if docking_metadata else "implemented_not_run"), "evidence": f"{docking_metadata.get('outputs', {}).get('compound_count', 0)} 个分子 {docking_metadata.get('outputs', {}).get('pose_count', 0)} 个姿势" if docking_metadata else "固定版本 输入准备 运行与SDF解析代码已实现"},
        {"step": 10, "name": "膜环境分子动力学", "status": "pending", "evidence": "需最终候选 膜体系参数与服务器"},
        {"step": 11, "name": "MM GBSA", "status": "pending", "evidence": "依赖合格 MD 轨迹"},
        {"step": 12, "name": "湿实验", "status": "external", "evidence": "不属于本仓库计算代码的完成范围"},
    ]

    return {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": purpose,
        "scientific_status": "computational_prediction_only",
        "active_pdb_ids": active_ids,
        "reference_only_pdb_ids": consensus_rules.get("reference_only_pdb_ids", []),
        "tier_definition": {str(key): value for key, value in sorted(tier_map.items())},
        "structures": structures,
        "total_consensus_regions": total_consensus,
        "alerts": alerts,
        "stages": stages,
        "drugclip": {
            "ranked_csv": str(ranked_csv.resolve()) if ranked_csv else "",
            "ranked_count": len(ranked_rows),
            "rows": ranked_rows,
            "input_run": str(input_run.resolve()) if input_run else "",
            "input_metadata": input_metadata,
            "selected_consensus": selected_consensus,
        },
        "docking": {
            "run": str(docking_run.resolve()) if docking_run else "",
            "metadata": docking_metadata,
            "best_poses": docking_best_rows,
        },
        "assessment": {
            "run": str(assessment_run.resolve()) if assessment_run else "",
            "metadata": assessment_metadata,
        },
        "stability": stability_metadata,
        "evidence_hashes": {
            "structures_config": sha256_file(root / "config" / "structures.json"),
            "consensus_rules": sha256_file(root / "config" / "consensus_rules.json"),
            "structure_manifest": sha256_file(manifest_path) if manifest_path.is_file() else "",
            "ranked_csv": sha256_file(ranked_csv) if ranked_csv and ranked_csv.is_file() else "",
            "input_run": sha256_file(input_run) if input_run and input_run.is_file() else "",
            "docking_run": sha256_file(docking_run) if docking_run and docking_run.is_file() else "",
            "assessment_run": sha256_file(assessment_run) if assessment_run and assessment_run.is_file() else "",
            "stability_json": sha256_file(stability_json) if stability_json and stability_json.is_file() else "",
            **{
                f"{pdb_id}_consensus_run": sha256_file(
                    root / "results" / "consensus" / f"{pdb_id}_consensus_pockets.run.json"
                )
                for pdb_id in active_ids
                if (root / "results" / "consensus" / f"{pdb_id}_consensus_pockets.run.json").is_file()
            },
        },
    }


def _status_text(status: str) -> str:
    return {
        "complete": "已完成",
        "validated": "已实跑验证",
        "complete_run": "已完成一次运行",
        "technical_validation_only": "仅工程验证",
        "implemented_not_run": "代码已实现 未附本次运行",
        "policy_only": "规则已固化 未运行",
        "pending": "待完成",
        "incomplete": "不完整",
        "external": "外部实验环节",
    }.get(status, status)


def render_html(snapshot: dict[str, Any]) -> str:
    esc = lambda value: html.escape(str(value), quote=True)
    structure_rows = []
    for item in snapshot["structures"]:
        counts = item["tool_counts"]
        consensus = item["consensus"]
        structure_rows.append(
            "<tr>"
            f"<td><strong>{esc(item['pdb_id'])}</strong></td>"
            f"<td>{esc(item['mutation_label'] or 'none')}</td>"
            f"<td>{esc(item['resolution_angstrom'])}</td>"
            f"<td>{counts['p2rank']}</td><td>{counts['dogsite3']}</td><td>{counts['fpocket']}</td>"
            f"<td>{consensus['count']}</td>"
            f"<td>{consensus['tier_counts']['T1']}</td><td>{consensus['tier_counts']['T2']}</td>"
            f"<td>{esc('; '.join(f'{key}={value}' for key, value in sorted(consensus['scope_counts'].items())))}</td>"
            "</tr>"
        )
    stage_rows = "".join(
        "<tr>"
        f"<td>{stage['step']}</td><td>{esc(stage['name'])}</td>"
        f"<td><span class='status status-{esc(stage['status'])}'>{esc(_status_text(stage['status']))}</span></td>"
        f"<td>{esc(stage['evidence'])}</td>"
        "</tr>"
        for stage in snapshot["stages"]
    )
    alert_rows = "".join(
        f"<li class='alert alert-{esc(alert['severity'])}'><strong>{esc(alert['title'])}</strong><span>{esc(alert['detail'])}</span></li>"
        for alert in snapshot["alerts"]
    ) or "<li class='alert'><strong>未发现结构级硬错误</strong><span>仍需进行生物学人工复核。</span></li>"

    ranked_rows = snapshot["drugclip"]["rows"][:20]
    scores = [float(row["drugclip_score"]) for row in ranked_rows if row.get("drugclip_score")]
    minimum = min(scores) if scores else 0.0
    maximum = max(scores) if scores else 0.0
    span = maximum - minimum
    result_rows = []
    for row in ranked_rows:
        score = float(row.get("drugclip_score", 0.0))
        width = 100.0 if span == 0 else 12.0 + 88.0 * (score - minimum) / span
        source_id = row.get("source_id", "")
        source_url = row.get("source_url", "")
        source_cell = esc(source_id or row.get("source", ""))
        if source_url:
            source_cell = f"<a href='{esc(source_url)}'>{source_cell}</a>"
        result_rows.append(
            "<tr>"
            f"<td>{esc(row.get('rank', ''))}</td>"
            f"<td><strong>{esc(row.get('compound_id', ''))}</strong></td>"
            f"<td>{source_cell}</td>"
            f"<td class='score'>{score:.6f}<div class='bar'><span style='width:{width:.1f}%'></span></div></td>"
            "</tr>"
        )
    result_section = (
        "<h2>DrugCLIP 排名结果</h2>"
        "<p>分数仅用于本次同一口袋与同一分子库内的相对排序，数值更高表示模型相似度更高；"
        "它不是结合自由能、药效、毒性或临床有效性的证据。</p>"
        "<div class='table-wrap'><table><thead><tr><th>排名</th><th>化合物</th><th>数据来源</th><th>DrugCLIP 分数</th></tr></thead>"
        f"<tbody>{''.join(result_rows)}</tbody></table></div>"
        if result_rows
        else "<h2>DrugCLIP 排名结果</h2><p>本报告未附带 DrugCLIP 排名文件。</p>"
    )
    selected = snapshot.get("drugclip", {}).get("selected_consensus", {})
    selected_section = ""
    if selected:
        selected_section = (
            "<h2>本次演示口袋</h2>"
            "<div class='table-wrap'><table><thead><tr><th>结构/口袋</th><th>层级</th><th>工具数</th><th>核心残基</th><th>最大中心距离 Å</th><th>最小 Jaccard 相似度</th><th>最大 Jaccard 距离</th></tr></thead>"
            "<tbody><tr>"
            f"<td><strong>{esc(selected.get('pdb_id', ''))}/{esc(selected.get('consensus_id', ''))}</strong></td>"
            f"<td>{esc(selected.get('tier', ''))}</td><td>{esc(selected.get('support_count', ''))}</td>"
            f"<td>{esc(selected.get('core_residue_count_2plus', ''))}</td>"
            f"<td>{esc(selected.get('max_member_center_distance', ''))}</td>"
            f"<td>{esc(selected.get('min_pairwise_jaccard_similarity', ''))}</td>"
            f"<td>{esc(selected.get('max_pairwise_jaccard_distance', ''))}</td>"
            "</tr></tbody></table></div>"
            "<p>Jaccard 相似度越接近 1 越相似；Jaccard 距离越接近 0 越相似。二者是同一信息的相反表达。</p>"
        )
    docking = snapshot.get("docking", {}).get("metadata", {})
    docking_section = ""
    if docking:
        outputs = docking.get("outputs", {})
        pose_rows = snapshot.get("docking", {}).get("best_poses", [])
        pose_table_rows = "".join(
            "<tr>"
            f"<td>{esc(row.get('docking_rank', ''))}</td>"
            f"<td><strong>{esc(row.get('compound_id', ''))}</strong></td>"
            f"<td>{esc(row.get('drugclip_rank', ''))}</td>"
            f"<td>{esc(row.get('CNNscore', ''))}</td>"
            f"<td>{esc(row.get('minimizedAffinity', ''))}</td>"
            "</tr>"
            for row in pose_rows
        )
        docking_section = (
            "<h2>GNINA docking 工程验证</h2>"
            f"<p>对 <strong>{esc(docking.get('pdb_id', ''))}/{esc(docking.get('consensus_id', ''))}</strong> "
            f"完成 {esc(outputs.get('compound_count', 0))} 个分子、{esc(outputs.get('pose_count', 0))} 个姿势的可审计运行。"
            "这些是刚性受体条件下的计算姿势与评分，不是结合或药效实验证据。</p>"
            + (
                "<div class='table-wrap'><table><thead><tr><th>Docking 排名</th><th>化合物</th><th>DrugCLIP 排名</th><th>CNNscore</th><th>minimizedAffinity</th></tr></thead>"
                f"<tbody>{pose_table_rows}</tbody></table></div>"
                if pose_table_rows else ""
            )
        )
    stability = snapshot.get("stability", {})
    stability_section = ""
    if stability:
        stability_section = (
            "<h2>数值精度稳定性</h2>"
            f"<p>同一批 {esc(stability.get('compound_count', 0))} 个分子的 FP16/FP32 排名 Spearman 相关系数为 "
            f"<strong>{float(stability.get('spearman_rank_correlation', 0)):.6f}</strong>，"
            f"Top {esc(stability.get('top_k', 0))} 重合 {esc(stability.get('top_k_overlap_count', 0))}/{esc(stability.get('top_k', 0))}。"
            "这只说明数值精度下排序稳定，不是生物学验证。</p>"
        )
    purpose_text = "工程验证数据 不可作为候选药结论" if snapshot["purpose"] != "formal_screening" else "指定分子库计算结果 仍需后续验证"
    return f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>PIEZO1 虚拟筛选证据报告</title>
<style>
:root{{--ink:#172033;--muted:#65728a;--line:#dce3ee;--blue:#1f5eff;--blue-soft:#eef4ff;--green:#087a55;--amber:#a05a00;--red:#b42318;--bg:#f5f7fb}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--bg);color:var(--ink);font-family:"Segoe UI","Microsoft YaHei",sans-serif;line-height:1.65}}
.page{{max-width:1180px;margin:0 auto;padding:42px 28px 70px}} header{{background:#fff;border:1px solid var(--line);border-radius:18px;padding:30px 34px;box-shadow:0 10px 30px rgba(26,45,85,.06)}}
h1{{margin:0 0 10px;font-size:32px;letter-spacing:-.4px}} h2{{margin:34px 0 12px;font-size:22px}} p{{margin:8px 0;color:#344057}}
.eyebrow{{color:var(--blue);font-weight:700;font-size:13px;letter-spacing:.08em;text-transform:uppercase}} .badge{{display:inline-block;margin-top:12px;padding:6px 11px;border-radius:999px;background:#fff4e5;color:#8a4b00;font-weight:700;font-size:13px}}
.cards{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:14px;margin-top:22px}} .card{{background:#fff;border:1px solid var(--line);border-radius:14px;padding:18px}} .metric{{font-size:29px;font-weight:750}} .label{{color:var(--muted);font-size:13px}}
.table-wrap{{overflow:auto;background:#fff;border:1px solid var(--line);border-radius:14px}} table{{width:100%;border-collapse:collapse;min-width:760px}} th,td{{padding:12px 14px;border-bottom:1px solid var(--line);text-align:left;vertical-align:middle}} th{{background:#edf2fb;font-size:13px}} tbody tr:last-child td{{border-bottom:0}}
.status{{display:inline-block;padding:4px 9px;border-radius:999px;font-size:12px;font-weight:700;white-space:nowrap;background:#eef1f6;color:#46536a}} .status-complete,.status-validated,.status-complete_run{{background:#e8f7f0;color:var(--green)}} .status-technical_validation_only,.status-policy_only,.status-implemented_not_run{{background:#fff4e5;color:var(--amber)}} .status-pending,.status-incomplete{{background:#fff0ee;color:var(--red)}}
.alerts{{list-style:none;padding:0;display:grid;gap:10px}} .alert{{background:#fff;border:1px solid var(--line);border-left:5px solid #98a2b3;border-radius:10px;padding:13px 16px}} .alert strong,.alert span{{display:block}} .alert span{{color:#4b5870;margin-top:3px}} .alert-high,.alert-critical{{border-left-color:var(--red)}}
.score{{font-variant-numeric:tabular-nums;min-width:220px}} .bar{{height:5px;background:#e7ecf5;border-radius:6px;margin-top:5px;overflow:hidden}} .bar span{{display:block;height:100%;background:linear-gradient(90deg,#72a0ff,var(--blue));border-radius:6px}} a{{color:var(--blue);text-decoration:none}} .foot{{margin-top:30px;color:var(--muted);font-size:12px}}
@media(max-width:780px){{.cards{{grid-template-columns:repeat(2,1fr)}} .page{{padding:20px 14px}} header{{padding:24px}}}}
</style></head><body><main class="page">
<header><div class="eyebrow">Reproducible evidence report</div><h1>PIEZO1 虚拟筛选证据报告</h1>
<p>报告把结构、三工具共识、DrugCLIP 排名和后续门禁分开呈现。所有分数均为计算预测，不能替代 docking、MD、毒理、湿实验或临床验证。</p><span class="badge">{esc(purpose_text)}</span></header>
<section class="cards"><div class="card"><div class="metric">{len(snapshot['active_pdb_ids'])}</div><div class="label">正式独立结构</div></div><div class="card"><div class="metric">3</div><div class="label">口袋预测工具</div></div><div class="card"><div class="metric">{snapshot['total_consensus_regions']}</div><div class="label">同结构共识区域</div></div><div class="card"><div class="metric">{snapshot['drugclip']['ranked_count']}</div><div class="label">本报告 DrugCLIP 分子数</div></div></section>
<h2>十二步方案执行状态</h2><div class="table-wrap"><table><thead><tr><th>Step</th><th>环节</th><th>状态</th><th>可核查证据</th></tr></thead><tbody>{stage_rows}</tbody></table></div>
<h2>四个结构的口袋结果</h2><p>T1 固定表示恰好 2 个工具支持，T2 固定表示 3 个工具支持。不同 PDB 的口袋没有合并。</p>
<div class="table-wrap"><table><thead><tr><th>PDB</th><th>突变</th><th>分辨率 Å</th><th>P2Rank</th><th>DoGSite3</th><th>fpocket</th><th>共识</th><th>T1</th><th>T2</th><th>口袋范围分布</th></tr></thead><tbody>{''.join(structure_rows)}</tbody></table></div>
<h2>必须保留的科学边界</h2><ul class="alerts">{alert_rows}</ul>{selected_section}{result_section}{docking_section}{stability_section}
<p class="foot">生成时间 UTC {esc(snapshot['generated_at_utc'])} · scientific_status={esc(snapshot['scientific_status'])}</p>
</main></body></html>"""


def write_report_bundle(output_dir: Path, snapshot: dict[str, Any]) -> dict[str, Path]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty report directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "summary.json"
    html_path = output_dir / "index.html"
    structure_path = output_dir / "structure_summary.csv"
    write_json(summary_path, snapshot)
    html_path.write_text(render_html(snapshot), encoding="utf-8")
    with structure_path.open("w", encoding="utf-8-sig", newline="") as handle:
        fields = [
            "pdb_id",
            "mutation_label",
            "resolution_angstrom",
            "p2rank_pockets",
            "dogsite3_pockets",
            "fpocket_pockets",
            "consensus_regions",
            "tier1_support_2",
            "tier2_support_3",
            "mutation_coordinate_status",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for item in snapshot["structures"]:
            writer.writerow(
                {
                    "pdb_id": item["pdb_id"],
                    "mutation_label": item["mutation_label"],
                    "resolution_angstrom": item["resolution_angstrom"],
                    "p2rank_pockets": item["tool_counts"]["p2rank"],
                    "dogsite3_pockets": item["tool_counts"]["dogsite3"],
                    "fpocket_pockets": item["tool_counts"]["fpocket"],
                    "consensus_regions": item["consensus"]["count"],
                    "tier1_support_2": item["consensus"]["tier_counts"]["T1"],
                    "tier2_support_3": item["consensus"]["tier_counts"]["T2"],
                    "mutation_coordinate_status": item["mutation_coordinate_status"],
                }
            )
    manifest_path = output_dir / "report_manifest.json"
    write_json(
        manifest_path,
        {
            "schema_version": 1,
            "files": {
                path.name: {"sha256": sha256_file(path), "bytes": path.stat().st_size}
                for path in (html_path, summary_path, structure_path)
            },
        },
    )
    return {
        "html": html_path,
        "summary": summary_path,
        "structures": structure_path,
        "manifest": manifest_path,
    }
