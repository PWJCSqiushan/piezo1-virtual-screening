from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable


TRUE_VALUES = {"1", "true", "yes", "y", "pass", "passed"}
FALSE_VALUES = {"0", "false", "no", "n", "fail", "failed"}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def normalized_key(row: dict[str, str]) -> str:
    compound_id = (row.get("compound_id") or "").strip()
    smiles = (row.get("canonical_smiles") or row.get("smiles") or "").strip()
    if compound_id:
        return f"id:{compound_id}"
    if smiles:
        return f"smiles:{smiles}"
    raise ValueError("Every assessment row needs compound_id or canonical_smiles")


def index_tool_rows(rows: list[dict[str, str]], tool: str) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        key = normalized_key(row)
        if key in indexed:
            raise ValueError(f"Duplicate {tool} assessment key: {key}")
        indexed[key] = row
    return indexed


def prefixed_values(row: dict[str, str] | None, prefix: str) -> dict[str, str]:
    if not row:
        return {}
    return {
        f"{prefix}__{key.strip()}": value
        for key, value in row.items()
        if key and key.strip() not in {"compound_id", "canonical_smiles", "smiles"}
    }


def assessment_branch(value: str | None) -> tuple[str, str]:
    normalized = (value or "").strip().lower()
    if normalized in TRUE_VALUES:
        return "oral_priority", "oral heuristic passed; expert review still required"
    if normalized in FALSE_VALUES:
        return (
            "inhalation_expert_review",
            "oral heuristic failed; retained for COPD inhalation and formulation expert review",
        )
    return "general_expert_review", "oral heuristic missing, uncertain or unrecognized"
