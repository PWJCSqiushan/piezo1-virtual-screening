from __future__ import annotations

import csv
import math
from pathlib import Path

from piezo_vs.drugclip_io import PocketAtom


CANDIDATE_LIST_REQUIRED_FIELDS = (
    "compound_id",
    "site_direction",
    "pdb_id",
    "consensus_id",
)
CANDIDATE_LIST_REVIEW_FIELDS = (
    "review_status",
    "approved",
    "evaluation_complete",
)
APPROVED_REVIEW_STATUSES = {"approved", "accepted", "manual_approved"}
TRUE_VALUES = {"1", "true", "yes", "y", "on"}
FALSE_VALUES = {"0", "false", "no", "n", "off"}

# The default policy removes every HETATM record.  Retaining waters/ions is a
# separate, reviewable choice and is recorded in the input manifest.
WATER_RESNAMES = {"HOH", "WAT", "DOD", "H2O"}
COMMON_ION_RESNAMES = {
    "CA",
    "CL",
    "CS",
    "FE",
    "K",
    "MG",
    "MN",
    "NA",
    "RB",
    "ZN",
}
HETATM_KEEP_RESNAMES = WATER_RESNAMES | COMMON_ION_RESNAMES
HETATM_POLICIES = {"remove_all", "keep_water_ions", "keep_all"}


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


def parse_bool(value: object, *, field: str) -> bool | None:
    """Parse a CSV boolean without silently treating an unknown value as false."""

    text = str(value or "").strip().lower()
    if not text:
        return None
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    raise ValueError(f"Invalid boolean for {field}: {value!r}")


def _read_csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        raw_fields = reader.fieldnames or []
        fields: list[str] = []
        field_by_lower: dict[str, str] = {}
        for raw_field in raw_fields:
            field = (raw_field or "").strip()
            key = field.lower()
            if not key:
                raise ValueError(f"CSV {path} contains an empty header")
            if key in field_by_lower:
                raise ValueError(f"CSV {path} contains duplicate header {field!r}")
            fields.append(key)
            field_by_lower[key] = field
        rows = []
        for row in reader:
            normalized = {
                field: (row.get(original) or "").strip()
                for field, original in field_by_lower.items()
            }
            # DictReader puts surplus cells under a None key.  A shifted CSV
            # must be rejected because row identity and SMILES would be unsafe.
            if None in row and row[None]:
                raise ValueError(f"CSV {path} has more cells than headers at row {reader.line_num}")
            rows.append(normalized)
    return fields, rows


def _canonicalize_candidate_fields(row: dict[str, str]) -> dict[str, str]:
    row = dict(row)
    aliases = {
        "source_rank": "rank",
        "source_drugclip_score": "drugclip_score",
        "smiles": "canonical_smiles",
        "smi": "canonical_smiles",
        "molecule_id": "compound_id",
        "id": "compound_id",
        "review": "review_status",
        "reviewer_status": "review_status",
        "evaluation_status": "evaluation_complete",
    }
    for source, target in aliases.items():
        if not row.get(target) and row.get(source):
            row[target] = row[source]
    for field in CANDIDATE_LIST_REQUIRED_FIELDS:
        row[field] = (row.get(field) or "").strip()
    review_status = row.get("review_status", "").strip().lower()
    row["review_status"] = review_status
    approved = parse_bool(row.get("approved"), field="approved")
    evaluation_complete = parse_bool(
        row.get("evaluation_complete"), field="evaluation_complete"
    )
    if approved is None and review_status in APPROVED_REVIEW_STATUSES:
        approved = True
    row["approved"] = "true" if approved is True else "false" if approved is False else ""
    row["evaluation_complete"] = (
        "true"
        if evaluation_complete is True
        else "false"
        if evaluation_complete is False
        else ""
    )
    return row


def read_candidate_list(
    path: Path,
    pdb_id: str,
    consensus_id: str,
    *,
    purpose: str = "technical_validation",
    ranked_rows: list[dict[str, str]] | None = None,
) -> list[dict[str, str]]:
    """Read and validate the explicit GNINA candidate-list contract.

    Required columns are ``compound_id``, ``site_direction``, ``pdb_id`` and
    ``consensus_id``.  ``canonical_smiles`` may be supplied by the selector or
    joined from a matching ranked CSV.  Formal runs require every row to have
    both ``approved=true`` and ``evaluation_complete=true``.
    """

    fields, raw_rows = _read_csv_rows(path)
    missing = [field for field in CANDIDATE_LIST_REQUIRED_FIELDS if field not in fields]
    if missing:
        raise ValueError(
            f"Candidate list {path} is missing required columns: {', '.join(missing)}"
        )
    if not raw_rows:
        raise ValueError(f"Candidate list {path} is empty")
    expected_pdb = pdb_id.upper()
    expected_consensus = consensus_id.upper()
    ranked_by_id = {
        row.get("compound_id", "").strip(): row for row in (ranked_rows or [])
    }
    candidates: list[dict[str, str]] = []
    seen: set[str] = set()
    for row_number, raw_row in enumerate(raw_rows, 2):
        row = _canonicalize_candidate_fields(raw_row)
        compound_id = row["compound_id"]
        if not compound_id:
            raise ValueError(f"Candidate list row {row_number} has an empty compound_id")
        if compound_id in seen:
            raise ValueError(f"Candidate list has duplicate compound_id {compound_id!r}")
        seen.add(compound_id)
        if not row["site_direction"]:
            raise ValueError(f"Candidate list row {row_number} has an empty site_direction")
        if row["pdb_id"].upper() != expected_pdb:
            raise ValueError(
                f"Candidate list row {row_number} has pdb_id={row['pdb_id']!r}; "
                f"expected {expected_pdb!r}"
            )
        if row["consensus_id"].upper() != expected_consensus:
            raise ValueError(
                f"Candidate list row {row_number} has consensus_id={row['consensus_id']!r}; "
                f"expected {expected_consensus!r}"
            )
        row["pdb_id"] = expected_pdb
        row["consensus_id"] = expected_consensus
        ranked = ranked_by_id.get(compound_id)
        candidate_smiles = row.get("canonical_smiles", "").strip()
        ranked_smiles = (ranked or {}).get("canonical_smiles", "").strip()
        if candidate_smiles and ranked_smiles and candidate_smiles != ranked_smiles:
            raise ValueError(
                f"Candidate {compound_id!r} has a canonical_smiles conflict with ranked CSV"
            )
        if not candidate_smiles:
            candidate_smiles = ranked_smiles
            row["canonical_smiles"] = candidate_smiles
        if not candidate_smiles:
            raise ValueError(
                f"Candidate {compound_id!r} has no canonical_smiles and is absent from ranked CSV"
            )
        if ranked:
            # Fill provenance fields from ranked output while preserving review
            # decisions and explicit site attribution from the selector.
            for field, value in ranked.items():
                if not row.get(field):
                    row[field] = value
        if purpose == "formal_screening":
            direction = "C006_C007" if expected_consensus in {"C006", "C007"} else expected_consensus
            if row["site_direction"] != direction:
                raise ValueError("Candidate direction does not match the docking site")
            approved = row.get("approved") == "true"
            evaluation_complete = row.get("evaluation_complete") == "true"
            review_ok = row.get("review_status", "") in APPROVED_REVIEW_STATUSES
            if not approved or not evaluation_complete or not review_ok:
                raise ValueError(
                    "Formal GNINA requires every candidate to be fully evaluated and manually "
                    f"approved; row {row_number} {compound_id!r} has "
                    f"review_status={row.get('review_status')!r}, "
                    f"approved={row.get('approved')!r}, "
                    f"evaluation_complete={row.get('evaluation_complete')!r}"
                )
        candidates.append(row)
    return candidates


def write_receptor_pdb(
    source: Path,
    output: Path,
    allowed_chains: set[str],
    *,
    heteroatom_policy: str = "remove_all",
) -> int:
    """Write the configured protein environment and return kept ATOM count.

    ``heteroatom_policy`` is independent from pocket extraction.  The legacy
    default removes all HETATM records; alternate policies are recorded in the
    input manifest by the caller and require review for formal runs.
    """

    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    if heteroatom_policy not in HETATM_POLICIES:
        raise ValueError(
            f"Unsupported heteroatom_policy={heteroatom_policy!r}; "
            f"choose one of {sorted(HETATM_POLICIES)}"
        )
    lines = ["HEADER    PIEZO1 RIGID RECEPTOR FOR COMPUTATIONAL DOCKING"]
    count = 0
    kept_hetatm = 0
    with source.open(encoding="ascii", errors="replace") as handle:
        for raw in handle:
            record = raw[:6]
            if record not in {"ATOM  ", "HETATM"}:
                continue
            chain = raw[21].strip()
            if record == "ATOM  " and chain not in allowed_chains:
                continue
            if raw[16].strip() not in {"", "A"}:
                continue
            if record == "ATOM  ":
                lines.append(raw.rstrip("\r\n"))
                count += 1
                continue
            resname = raw[17:20].strip().upper()
            keep = heteroatom_policy == "keep_all" or (
                heteroatom_policy == "keep_water_ions" and resname in HETATM_KEEP_RESNAMES
            )
            if keep:
                lines.append(raw.rstrip("\r\n"))
                kept_hetatm += 1
    if count == 0:
        raise ValueError(f"No receptor ATOM records found for chains {sorted(allowed_chains)}")
    lines.extend(
        [
            f"REMARK 900 HETATM POLICY {heteroatom_policy}",
            f"REMARK 900 HETATM KEPT {kept_hetatm}",
            "REMARK 900 PRIMARY ALTLOC BLANK OR A",
            "REMARK 900 PROTONATION AND MEMBRANE ENVIRONMENT REQUIRE EXPERT REVIEW",
            "END",
        ]
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="ascii", errors="replace")
    return count


def receptor_policy_summary(
    source: Path,
    allowed_chains: set[str],
    *,
    heteroatom_policy: str = "remove_all",
) -> dict[str, object]:
    """Summarize protein and HETATM decisions without writing a receptor."""

    if heteroatom_policy not in HETATM_POLICIES:
        raise ValueError(
            f"Unsupported heteroatom_policy={heteroatom_policy!r}; "
            f"choose one of {sorted(HETATM_POLICIES)}"
        )
    atom_count = 0
    selected_hetatm = 0
    kept_hetatm = 0
    removed_hetatm = 0
    kept_resnames: dict[str, int] = {}
    with source.open(encoding="ascii", errors="replace") as handle:
        for raw in handle:
            record = raw[:6]
            if record not in {"ATOM  ", "HETATM"}:
                continue
            if (record == "ATOM  " and raw[21].strip() not in allowed_chains) or raw[16].strip() not in {"", "A"}:
                continue
            if record == "ATOM  ":
                atom_count += 1
                continue
            selected_hetatm += 1
            resname = raw[17:20].strip().upper()
            keep = heteroatom_policy == "keep_all" or (
                heteroatom_policy == "keep_water_ions" and resname in HETATM_KEEP_RESNAMES
            )
            if keep:
                kept_hetatm += 1
                kept_resnames[resname] = kept_resnames.get(resname, 0) + 1
            else:
                removed_hetatm += 1
    return {
        "protein_atom_count": atom_count,
        "selected_hetatm_count": selected_hetatm,
        "kept_hetatm_count": kept_hetatm,
        "removed_hetatm_count": removed_hetatm,
        "kept_hetatm_resnames": dict(sorted(kept_resnames.items())),
        "heteroatom_policy": heteroatom_policy,
        "water_resnames": sorted(WATER_RESNAMES),
        "ion_resnames": sorted(COMMON_ION_RESNAMES),
    }


def read_ranked_compounds(
    path: Path,
    pdb_id: str,
    consensus_id: str,
    top_n: int | None,
) -> list[dict[str, str]]:
    fields, rows = _read_csv_rows(path)
    if not rows:
        raise ValueError("Ranked compound CSV is empty")
    if top_n is not None and top_n < 1:
        raise ValueError("top_n must be positive")
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
    return selected if top_n is None else selected[:top_n]


def numeric_property(properties: dict[str, object], name: str) -> float | None:
    value = properties.get(name, "")
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None
