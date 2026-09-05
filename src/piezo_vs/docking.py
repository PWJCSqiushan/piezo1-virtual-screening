from __future__ import annotations

import csv
from pathlib import Path

from piezo_vs.drugclip_io import PocketAtom


def docking_box_from_atoms(
    atoms: list[PocketAtom],
    center: tuple[float, float, float],
    *,
    padding: float = 6.0,
    minimum_size: float = 18.0,
    maximum_size: float = 40.0,
) -> dict[str, object]:
    if not atoms:
        raise ValueError("Cannot build a docking box without pocket atoms")
    raw_sizes = []
    sizes = []
    for axis in range(3):
        radius = max(abs(atom.coordinate[axis] - center[axis]) for atom in atoms)
        raw = 2.0 * radius + 2.0 * padding
        raw_sizes.append(raw)
        sizes.append(min(max(raw, minimum_size), maximum_size))
    return {
        "center_x": round(center[0], 3),
        "center_y": round(center[1], 3),
        "center_z": round(center[2], 3),
        "size_x": round(sizes[0], 3),
        "size_y": round(sizes[1], 3),
        "size_z": round(sizes[2], 3),
        "raw_size_x": round(raw_sizes[0], 3),
        "raw_size_y": round(raw_sizes[1], 3),
        "raw_size_z": round(raw_sizes[2], 3),
        "padding": padding,
        "minimum_size": minimum_size,
        "maximum_size": maximum_size,
        "clipped_to_maximum": any(raw > maximum_size for raw in raw_sizes),
    }


def write_receptor_pdb(source: Path, output: Path, allowed_chains: set[str]) -> int:
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    lines = ["HEADER    PIEZO1 RIGID RECEPTOR FOR COMPUTATIONAL DOCKING"]
    count = 0
    with source.open(encoding="ascii", errors="replace") as handle:
        for raw in handle:
            if not raw.startswith("ATOM  "):
                continue
            if raw[21].strip() not in allowed_chains:
                continue
            if raw[16].strip() not in {"", "A"}:
                continue
            lines.append(raw.rstrip("\r\n"))
            count += 1
    if count == 0:
        raise ValueError(f"No receptor ATOM records found for chains {sorted(allowed_chains)}")
    lines.extend(
        [
            "REMARK 900 HETATM RECORDS REMOVED AND PRIMARY ALTLOC RETAINED",
            "REMARK 900 PROTONATION AND MEMBRANE ENVIRONMENT REQUIRE EXPERT REVIEW",
            "END",
        ]
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="ascii", errors="replace")
    return count


def read_ranked_compounds(path: Path, pdb_id: str, consensus_id: str, top_n: int) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("Ranked compound CSV is empty")
    selected = []
    seen: set[str] = set()
    for row in rows:
        if (row.get("pdb_id") or "").upper() != pdb_id.upper():
            raise ValueError("Ranked CSV contains a different PDB ID")
        if (row.get("consensus_id") or "").upper() != consensus_id.upper():
            raise ValueError("Ranked CSV contains a different consensus ID")
        compound_id = row.get("compound_id", "")
        if not compound_id or compound_id in seen:
            raise ValueError(f"Missing or duplicate compound_id: {compound_id!r}")
        seen.add(compound_id)
        selected.append(row)
        if len(selected) == top_n:
            break
    return selected


def numeric_property(properties: dict[str, str], name: str) -> float | None:
    value = properties.get(name, "")
    try:
        return float(value)
    except (TypeError, ValueError):
        return None
