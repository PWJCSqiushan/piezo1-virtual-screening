from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.io_utils import read_json, sha256_file, write_json
from piezo_vs.pockets import (
    build_consensus_components,
    classify_pocket_scope,
    group_overlapping_consensus_components,
    parse_dogsite,
    parse_fpocket,
    parse_p2rank,
    summarize_component,
    validate_fpocket_output_name,
)


def latest_dogsite_downloads(pdb_id: str) -> Path:
    base = ROOT / "runs" / "dogsite3" / pdb_id
    latest = read_json(base / "latest_job.json")
    return base / latest["job_id"] / "downloads"


def latest_fpocket_output(pdb_id: str) -> Path:
    latest = read_json(ROOT / "runs" / "fpocket" / pdb_id / "latest_run.json")
    return ROOT / latest["output_dir"]


def parse_limit(value: str) -> int | None:
    if value.lower() == "all":
        return None
    limit = int(value)
    if limit <= 0:
        raise argparse.ArgumentTypeError("top-per-tool must be a positive integer or 'all'")
    return limit


def validate_fpocket_output(pdb_id: str, output_dir: Path) -> Path:
    output_dir = output_dir.resolve()
    try:
        validate_fpocket_output_name(pdb_id, output_dir)
    except ValueError as exc:
        raise SystemExit(
            f"{exc}. This check prevents results from another PDB being mislabeled."
        ) from exc
    target_pdb = ROOT / "data" / "raw" / "structures" / f"{pdb_id}.pdb"
    run_metadata = output_dir.parent / "run.json"
    if run_metadata.exists():
        metadata = read_json(run_metadata)
        if metadata.get("pdb_id") != pdb_id:
            raise SystemExit(f"fpocket run.json belongs to {metadata.get('pdb_id')}, not {pdb_id}")
        expected_hash = metadata.get("input_sha256")
        if expected_hash and expected_hash != sha256_file(target_pdb):
            raise SystemExit("fpocket input SHA256 does not match the configured target PDB")
    return output_dir


def public_source_path(source_file: str) -> str:
    path = Path(source_file).resolve()
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return f"external/{path.name}"


def main() -> int:
    parser = argparse.ArgumentParser(description="Build same-structure multi-tool pocket consensus.")
    parser.add_argument("--pdb-id", required=True)
    parser.add_argument("--fpocket-output", type=Path)
    parser.add_argument("--top-per-tool", default="all", help="Positive integer or 'all' (default).")
    parser.add_argument("--max-center-distance", type=float, default=12.0)
    parser.add_argument("--min-shared-residues", type=int, default=3)
    parser.add_argument("--allow-partial", action="store_true", help="Write preliminary 2-tool regions if fpocket is unavailable.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "results" / "consensus",
        help="Output directory; defaults to results/consensus.",
    )
    args = parser.parse_args()

    pdb_id = args.pdb_id.upper()
    limit = parse_limit(args.top_per_tool)
    rules = read_json(ROOT / "config" / "consensus_rules.json")
    if rules.get("forbid_cross_conformation_consensus") is not True:
        raise RuntimeError("Configuration must explicitly forbid cross-conformation consensus")
    if pdb_id not in rules["active_pdb_ids"]:
        raise SystemExit(f"{pdb_id} is not an active independently analyzed structure")

    p2rank_path = ROOT / "runs" / "p2rank" / pdb_id / "output" / f"{pdb_id}.pdb_predictions.csv"
    structures = read_json(ROOT / "config" / "structures.json")
    structure = next(item for item in structures["structures"] if item["pdb_id"] == pdb_id)
    pockets = parse_p2rank(p2rank_path, limit)
    pockets.extend(parse_dogsite(latest_dogsite_downloads(pdb_id), pdb_id, limit))
    fpocket_output = args.fpocket_output
    if fpocket_output is None:
        try:
            fpocket_output = latest_fpocket_output(pdb_id)
        except FileNotFoundError:
            if not args.allow_partial:
                raise SystemExit(
                    f"fpocket output is missing for {pdb_id}. Run scripts/run_fpocket.py first, "
                    "pass --fpocket-output, or use --allow-partial for a non-final precheck."
                )
    if fpocket_output is not None:
        fpocket_output = validate_fpocket_output(pdb_id, fpocket_output)
        pockets.extend(parse_fpocket(fpocket_output, limit))

    pockets = [replace(pocket, source_file=public_source_path(pocket.source_file)) for pocket in pockets]

    available_tools = sorted({pocket.tool for pocket in pockets})
    required_tools = sorted(rules["required_tools"])
    is_final = available_tools == required_tools
    if not is_final and not args.allow_partial:
        raise SystemExit(f"Required tools={required_tools}; available tools={available_tools}")

    tier_by_support = {
        int(definition["support_count"]): tier for tier, definition in rules["tiers"].items()
    }
    hypotheses = build_consensus_components(
        pockets,
        max_center_distance=args.max_center_distance,
        min_shared_residues=args.min_shared_residues,
    )
    groups = group_overlapping_consensus_components(
        hypotheses,
        max_center_distance=args.max_center_distance,
        min_shared_residues=args.min_shared_residues,
    )
    components = [representative for representative, _ in groups]
    pocket_memberships = Counter(
        (pocket.tool, pocket.pocket_id, pocket.source_file)
        for component in components
        for pocket in component
    )
    rows = []
    hypothesis_rows = []
    for index, (component, group_members) in enumerate(groups, 1):
        summary = summarize_component(component)
        support_count = int(summary["support_count"])
        scope, chains = classify_pocket_scope(
            component,
            set(structure["auth_piezo1_chains"]),
            set(structure["auth_mdfic_chains"]),
        )
        rows.append(
            {
                "consensus_id": f"C{index:03d}",
                "pdb_id": pdb_id,
                "classification_status": "final_three_tool_run" if is_final else "preliminary_partial_run",
                "tier": tier_by_support.get(support_count, "") if is_final else "",
                "pocket_scope": scope,
                "chains": chains,
                "source_pocket_reuse_warning": any(
                    pocket_memberships[(pocket.tool, pocket.pocket_id, pocket.source_file)] > 1
                    for pocket in component
                ),
                "max_source_pocket_region_memberships": max(
                    pocket_memberships[(pocket.tool, pocket.pocket_id, pocket.source_file)]
                    for pocket in component
                ),
                "overlap_group_size": len(group_members),
                **summary,
            }
        )
        for hypothesis_index, hypothesis in enumerate(group_members, 1):
            hypothesis_rows.append(
                {
                    "consensus_id": f"C{index:03d}",
                    "hypothesis_index": hypothesis_index,
                    "is_representative": hypothesis is component,
                    **summarize_component(hypothesis),
                }
            )

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{pdb_id}_consensus_pockets" if is_final else f"{pdb_id}_preliminary_pockets"
    csv_path = output_dir / f"{stem}.csv"
    if rows:
        with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    else:
        csv_path.write_text("consensus_id,pdb_id,classification_status,tier,support_count,tools\n", encoding="utf-8")
    write_json(output_dir / f"{stem}.json", rows)
    audit_csv_path = output_dir / f"{stem}.hypotheses.csv"
    if hypothesis_rows:
        with audit_csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(hypothesis_rows[0]))
            writer.writeheader()
            writer.writerows(hypothesis_rows)
    write_json(output_dir / f"{stem}.hypotheses.json", hypothesis_rows)
    metadata_path = output_dir / f"{stem}.run.json"
    metadata = {
        "schema_version": 1,
        "consensus_algorithm": "maximal_strict_cliques_grouped_by_fixed_representative_v2",
        "source_pockets_may_support_multiple_regions": True,
        "hypothesis_count": len(hypotheses),
        "region_group_count": len(groups),
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "pdb_id": pdb_id,
        "available_tools": available_tools,
        "required_tools": required_tools,
        "classification_status": "final_three_tool_run" if is_final else "preliminary_partial_run",
        "top_per_tool": "all" if limit is None else limit,
        "max_center_distance": args.max_center_distance,
        "min_shared_residues": args.min_shared_residues,
        "rules_sha256": sha256_file(ROOT / "config" / "consensus_rules.json"),
        "structures_sha256": sha256_file(ROOT / "config" / "structures.json"),
        "output_csv_sha256": sha256_file(csv_path),
        "hypotheses_csv_sha256": sha256_file(audit_csv_path) if audit_csv_path.is_file() else "",
        "region_count": len(rows),
    }
    write_json(metadata_path, metadata)
    print(
        f"{pdb_id}: tools={','.join(available_tools)}; status={'final' if is_final else 'preliminary'}; "
        f"consensus_regions={len(rows)}"
    )
    print(f"Wrote {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
