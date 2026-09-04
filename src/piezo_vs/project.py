from __future__ import annotations

from pathlib import Path


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def ensure_project_dirs(root: Path) -> None:
    for relative in (
        "data/raw/structures/api",
        "data/processed/structures",
        "data/processed/pockets",
        "results",
        "runs",
        "tools/cache",
    ):
        (root / relative).mkdir(parents=True, exist_ok=True)
