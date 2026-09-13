"""Step 8 candidate assessment and deterministic diversity selection.

This module is intentionally independent from the historical Step 8 import
helpers in :mod:`piezo_vs.assessment`.  The new candidate batch has one
canonical master table and two local DrugCLIP ranking directions.  Every
input row is retained in the assessment table; structure flags are review
signals and never constitute an automatic chemical rejection.

RDKit is an optional import at module import time because the public
repository's light-weight test environment does not ship the WSL DrugCLIP
runtime.  A real assessment requires RDKit, and rows are marked pending when
it is unavailable or a SMILES cannot be parsed.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


DIRECTIONS = ("C006_C007", "C016")
TOOL_NAMES = ("swissadme", "admetlab3", "protox3")
MASTER_REQUIRED_FIELDS = ("compound_id", "canonical_smiles")
RANK_REQUIRED_FIELDS = ("compound_id", "rank", "drugclip_score")
MORGAN_RADIUS = 2
MORGAN_NBITS = 2048
MORGAN_USE_CHIRALITY = True
MAX_SELECTION_POOL = 50
DEFAULT_SELECTION_COUNT = 10
MAX_SWISSADME_BATCH = 200

_TRUE_VALUES = {"1", "true", "yes", "y", "approved", "approve", "pass", "passed"}
_FALSE_VALUES = {"0", "false", "no", "n", "rejected", "reject", "fail", "failed"}
_FAILURE_VALUES = {
    "failed",
    "failure",
    "error",
    "errored",
    "invalid",
    "timeout",
    "timed_out",
    "not_run",
    "unavailable",
}

# The allowed set is deliberately conservative for the structure review
# flag.  It does not remove a compound from the pool.  H is normally implicit
# in a SMILES and is included for explicitly written hydrogen atoms.
_COMMON_ORGANIC_ELEMENTS = {"H", "B", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}
_METAL_ATOMIC_NUMBERS = (
    set(range(3, 5))
    | {11, 12, 13, 19, 20}
    | set(range(21, 33))
    | {37, 38}
    | set(range(39, 51))
    | {55, 56}
    | set(range(57, 84))
    | {87, 88}
    | set(range(89, 113))
)


class CandidateValidationError(ValueError):
    """Raised when a candidate, ranking, or review input is unsafe to use."""


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    """Read a UTF-8 CSV while preserving every source column as text."""

    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = [name.strip() for name in (reader.fieldnames or []) if name]
        if not fieldnames:
            raise CandidateValidationError(f"CSV has no header: {path}")
        missing_headers = [name for name in fieldnames if name not in (reader.fieldnames or [])]
        if missing_headers:
            raise CandidateValidationError(f"CSV header contains blank/invalid names: {path}")
        rows: list[dict[str, str]] = []
        for row in reader:
            rows.append({str(key): (value or "") for key, value in row.items() if key is not None})
        return rows


def write_csv_rows(
    path: Path, rows: Iterable[Mapping[str, object]], fields: Sequence[str] | None = None
) -> None:
    """Write deterministic UTF-8 CSV output, creating only the parent path."""

    materialized = [dict(row) for row in rows]
    if fields is None:
        discovered: list[str] = []
        for row in materialized:
            for field in row:
                if field not in discovered:
                    discovered.append(field)
        fields = discovered
    fields = list(fields)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in materialized:
            writer.writerow({field: _csv_value(row.get(field, "")) for field in fields})


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _csv_value(value: object) -> object:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return ""
        return format(value, ".12g")
    if isinstance(value, (list, tuple, set)):
        return ";".join(str(item) for item in value)
    if isinstance(value, dict):
        return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return value


def _text(value: object) -> str:
    return "" if value is None else str(value).strip()


def _normal_bool(value: object) -> str:
    normalized = _text(value).lower()
    if normalized in _TRUE_VALUES:
        return "true"
    if normalized in _FALSE_VALUES:
        return "false"
    return ""


def _safe_float(value: object, *, field: str, row_number: int) -> float:
    text = _text(value)
    if not text:
        raise CandidateValidationError(f"Missing {field} at row {row_number}")
    try:
        parsed = float(text)
    except ValueError as exc:
        raise CandidateValidationError(f"Invalid {field} at row {row_number}: {text!r}") from exc
    if not math.isfinite(parsed):
        raise CandidateValidationError(f"Non-finite {field} at row {row_number}: {text!r}")
    return parsed


def _safe_rank(value: object, *, row_number: int) -> int:
    text = _text(value)
    if not re.fullmatch(r"[1-9][0-9]*", text):
        raise CandidateValidationError(f"Invalid rank at row {row_number}: {text!r}")
    rank = int(text)
    if rank < 1:
        raise CandidateValidationError(f"Rank must be positive at row {row_number}")
    return rank


def validate_master_rows(rows: Iterable[Mapping[str, object]]) -> list[dict[str, str]]:
    """Validate the de-duplicated master table.

    The master table is intentionally strict: one ``compound_id`` maps to
    one canonical SMILES.  Source rows and source task identifiers belong in
    extra columns and are retained unchanged.
    """

    materialized = [dict(row) for row in rows]
    if not materialized:
        raise CandidateValidationError("Candidate master is empty")
    for field in MASTER_REQUIRED_FIELDS:
        if field not in materialized[0]:
            raise CandidateValidationError(f"Candidate master needs {field!r}")
    seen: dict[str, str] = {}
    seen_structures: set[str] = set()
    normalized: list[dict[str, str]] = []
    for row_number, raw in enumerate(materialized, 2):
        row = {str(key): _text(value) for key, value in raw.items()}
        compound_id = row.get("compound_id", "")
        smiles = row.get("canonical_smiles", "")
        if not compound_id or not smiles:
            raise CandidateValidationError(
                f"Master row {row_number} needs compound_id and canonical_smiles"
            )
        if compound_id in seen:
            if seen[compound_id] != smiles:
                raise CandidateValidationError(
                    f"Duplicate compound_id with conflicting canonical_smiles at row {row_number}: "
                    f"{compound_id!r}"
                )
            raise CandidateValidationError(f"Duplicate compound_id at row {row_number}: {compound_id!r}")
        if smiles in seen_structures:
            raise CandidateValidationError("Duplicate structure in canonical master")
        if row.get("structure_sha256") and row["structure_sha256"] != hashlib.sha256(smiles.encode()).hexdigest():
            raise CandidateValidationError("Master structure hash mismatch")
        seen_structures.add(smiles)
        seen[compound_id] = smiles
        row["compound_id"] = compound_id
        row["canonical_smiles"] = smiles
        normalized.append(row)
    return normalized


def _direction_from_row(row: Mapping[str, object], fallback: str | None) -> str | None:
    value = _text(row.get("direction") or row.get("site_direction") or fallback)
    if not value:
        return None
    value = value.upper().replace("-", "_")
    aliases = {
        "C006+C007": "C006_C007",
        "C006/C007": "C006_C007",
        "C006_C007": "C006_C007",
        "C016": "C016",
    }
    return aliases.get(value, value)


def validate_ranking_rows(
    rows: Iterable[Mapping[str, object]],
    direction: str,
    *,
    master_ids: set[str] | None = None,
) -> list[dict[str, object]]:
    """Validate one complete local ranking for a direction.

    Ranks must be unique and contiguous from one.  Scores must be finite;
    this rejects NaN and infinity before they can influence sorting or
    selection.  A ranking can cover a subset of the master, but it cannot
    introduce an unknown identifier.
    """

    direction = _direction_from_row({}, direction) or ""
    if direction not in DIRECTIONS:
        raise CandidateValidationError(f"Unsupported direction: {direction!r}")
    materialized = [dict(row) for row in rows]
    if not materialized:
        raise CandidateValidationError(f"Ranking for {direction} is empty")
    for field in RANK_REQUIRED_FIELDS:
        if field not in materialized[0]:
            raise CandidateValidationError(f"Ranking for {direction} needs {field!r}")
    seen_ids: set[str] = set()
    seen_ranks: set[int] = set()
    normalized: list[dict[str, object]] = []
    for row_number, raw in enumerate(materialized, 2):
        row = {str(key): _text(value) for key, value in raw.items()}
        compound_id = row.get("compound_id", "")
        if not compound_id:
            raise CandidateValidationError(f"Missing compound_id in {direction} row {row_number}")
        if compound_id in seen_ids:
            raise CandidateValidationError(
                f"Duplicate compound_id in {direction} ranking at row {row_number}: {compound_id!r}"
            )
        rank = _safe_rank(row.get("rank"), row_number=row_number)
        if rank in seen_ranks:
            raise CandidateValidationError(
                f"Duplicate rank in {direction} ranking at row {row_number}: {rank}"
            )
        score = _safe_float(row.get("drugclip_score"), field="drugclip_score", row_number=row_number)
        direction_in_row = _direction_from_row(row, None)
        if direction_in_row and direction_in_row != direction:
            raise CandidateValidationError(
                f"Ranking direction mismatch at row {row_number}: expected {direction}, "
                f"got {direction_in_row}"
            )
        if master_ids is not None and compound_id not in master_ids:
            raise CandidateValidationError(
                f"Ranking {direction} contains unknown compound_id: {compound_id!r}"
            )
        seen_ids.add(compound_id)
        seen_ranks.add(rank)
        row["compound_id"] = compound_id
        row["rank"] = rank
        row["drugclip_score"] = score
        row["direction"] = direction
        normalized.append(row)
    expected_ranks = set(range(1, len(normalized) + 1))
    if seen_ranks != expected_ranks:
        missing = sorted(expected_ranks - seen_ranks)
        extra = sorted(seen_ranks - expected_ranks)
        raise CandidateValidationError(
            f"Ranks for {direction} must be contiguous 1..{len(normalized)}; "
            f"missing={missing}, out_of_range={extra}"
        )
    ordered = sorted(normalized, key=lambda row: (int(row["rank"]), str(row["compound_id"])))
    if master_ids is not None and seen_ids != master_ids:
        raise CandidateValidationError("Local ranking must cover the complete candidate universe")
    expected = sorted(ordered, key=lambda row: (-float(row["drugclip_score"]), str(row["compound_id"])))
    if [r["compound_id"] for r in ordered] != [r["compound_id"] for r in expected]:
        raise CandidateValidationError("Ranks contradict scores or stable tie order")
    return ordered


def validate_rankings(
    ranking_rows_by_direction: Mapping[str, Iterable[Mapping[str, object]]],
    *,
    master_ids: set[str] | None = None,
) -> dict[str, list[dict[str, object]]]:
    normalized: dict[str, list[dict[str, object]]] = {}
    for direction, rows in ranking_rows_by_direction.items():
        normalized_direction = _direction_from_row({}, direction) or ""
        if normalized_direction in normalized:
            raise CandidateValidationError(f"Duplicate ranking direction: {normalized_direction}")
        normalized[normalized_direction] = validate_ranking_rows(
            rows, normalized_direction, master_ids=master_ids
        )
    return normalized


def _rdkit_modules() -> tuple[Any, Any, Any, Any] | None:
    """Return RDKit modules lazily, or None in the light-weight environment."""

    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors
    except ImportError:
        return None
    return Chem, RDLogger, (AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors), None


def rdkit_available() -> bool:
    return _rdkit_modules() is not None


def _pains_catalog(Chem: Any) -> Any | None:
    try:
        from rdkit.Chem import FilterCatalog

        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        return FilterCatalog.FilterCatalog(params)
    except (AttributeError, ImportError, RuntimeError):
        return None


def _stereo_is_undefined(Chem: Any, mol: Any) -> bool:
    try:
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        if any(label == "?" for _, label in centers):
            return True
    except (AttributeError, RuntimeError):
        pass
    try:
        potential = Chem.FindPotentialStereo(mol)
        for info in potential:
            specified = str(getattr(info, "specified", "")).lower()
            if "unspec" in specified:
                return True
    except (AttributeError, RuntimeError):
        pass
    return False


def _is_metal(atom: Any) -> bool:
    try:
        return int(atom.GetAtomicNum()) in _METAL_ATOMIC_NUMBERS
    except (AttributeError, TypeError, ValueError):
        return False


def morgan_fingerprint(
    smiles: str,
    *,
    radius: int = MORGAN_RADIUS,
    nbits: int = MORGAN_NBITS,
    use_chirality: bool = MORGAN_USE_CHIRALITY,
) -> str:
    """Return a deterministic Morgan bit string for a valid SMILES.

    A bit string rather than a pickled RDKit object keeps CSV provenance
    portable between Windows and WSL.  The fixed parameters are included in
    every assessment row and selection manifest.
    """

    modules = _rdkit_modules()
    if modules is None:
        raise RuntimeError("RDKit is required to generate Morgan fingerprints")
    Chem, _, rdkit_chem, _ = modules
    AllChem = rdkit_chem[0]
    mol = Chem.MolFromSmiles(_text(smiles))
    if mol is None:
        raise CandidateValidationError(f"Cannot generate fingerprint for invalid SMILES: {smiles!r}")
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=int(radius), nBits=int(nbits), useChirality=bool(use_chirality)
    )
    return fp.ToBitString()


def assess_smiles(smiles: str) -> dict[str, object]:
    """Calculate RDKit descriptors and non-destructive structure flags."""

    text = _text(smiles)
    empty = {
        "rdkit_status": "parse_failed",
        "rdkit_canonical_smiles": "",
        "molecular_weight": "",
        "mol_wt": "",
        "logp": "",
        "hbd": "",
        "hba": "",
        "rotatable_bonds": "",
        "tpsa": "",
        "heavy_atoms": "",
        "formal_charge": "",
        "ring_count": "",
        "fragment_count": "",
        "element_set": "",
        "has_salt": "",
        "has_metal": "",
        "has_isotope": "",
        "stereo_undefined": "",
        "pains_alert": "",
        "pains_alert_count": "",
        "structure_flags": "rdkit_parse_failed",
        "morgan_radius": MORGAN_RADIUS,
        "morgan_nbits": MORGAN_NBITS,
        "morgan_use_chirality": MORGAN_USE_CHIRALITY,
        "morgan_fingerprint": "",
    }
    modules = _rdkit_modules()
    if modules is None:
        unavailable = dict(empty)
        unavailable["rdkit_status"] = "rdkit_unavailable"
        unavailable["structure_flags"] = "rdkit_unavailable"
        return unavailable
    Chem, RDLogger, rdkit_chem, _ = modules
    AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors = rdkit_chem
    # Invalid external SMILES should be a data status, not a noisy console
    # error.  RDKit logging is restored by the interpreter at process end.
    try:
        RDLogger.DisableLog("rdApp.error")
        mol = Chem.MolFromSmiles(text)
    except (RuntimeError, ValueError):
        mol = None
    if mol is None:
        return empty

    atoms = list(mol.GetAtoms())
    elements = sorted({atom.GetSymbol() for atom in atoms})
    fragment_count = len(Chem.GetMolFrags(mol))
    has_salt = fragment_count > 1
    has_metal = any(_is_metal(atom) for atom in atoms)
    has_isotope = any(int(atom.GetIsotope()) != 0 for atom in atoms)
    stereo_undefined = _stereo_is_undefined(Chem, mol)
    catalog = _pains_catalog(Chem)
    pains_count: int | None = None
    if catalog is not None:
        try:
            pains_count = len(list(catalog.GetMatches(mol)))
        except (AttributeError, RuntimeError):
            pains_count = None

    flags: list[str] = []
    if has_salt:
        flags.append("multiple_fragments_or_salt")
    if has_metal:
        flags.append("metal_present")
    nonstandard = sorted(set(elements) - _COMMON_ORGANIC_ELEMENTS)
    if nonstandard:
        flags.append("nonstandard_element:" + ",".join(nonstandard))
    if has_isotope:
        flags.append("isotope_labeled")
    if stereo_undefined:
        flags.append("stereochemistry_unspecified")
    if pains_count is None:
        flags.append("pains_check_unavailable")
    elif pains_count:
        flags.append("pains_alert")

    molecular_weight = float(Descriptors.MolWt(mol))
    logp = float(Crippen.MolLogP(mol))
    hbd = int(Lipinski.NumHDonors(mol))
    hba = int(Lipinski.NumHAcceptors(mol))
    rotatable = int(Lipinski.NumRotatableBonds(mol))
    tpsa = float(rdMolDescriptors.CalcTPSA(mol))
    heavy_atoms = int(mol.GetNumHeavyAtoms())
    formal_charge = int(Chem.GetFormalCharge(mol))
    ring_count = int(rdMolDescriptors.CalcNumRings(mol))
    if heavy_atoms < 6: flags.append("very_small_molecule")
    if "C" not in elements: flags.append("no_carbon_atoms")
    if molecular_weight < 100: flags.append("low_molecular_weight")
    try:
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(
            mol,
            radius=MORGAN_RADIUS,
            nBits=MORGAN_NBITS,
            useChirality=MORGAN_USE_CHIRALITY,
        ).ToBitString()
    except (AttributeError, RuntimeError, ValueError):
        fingerprint = ""
        flags.append("morgan_fingerprint_failed")

    values: dict[str, object] = {
        "rdkit_status": "parsed",
        "rdkit_canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
        "lipinski_violation_count": sum((molecular_weight > 500, logp > 5, hbd > 5, hba > 10)),
        "lipinski_rule_note": "RDKit MW/logP/HBD/HBA approximation; review signal, no automatic rejection or route claim",
        "molecular_weight": round(molecular_weight, 6),
        "mol_wt": round(molecular_weight, 6),
        "logp": round(logp, 6),
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rotatable,
        "tpsa": round(tpsa, 6),
        "heavy_atoms": heavy_atoms,
        "formal_charge": formal_charge,
        "ring_count": ring_count,
        "fragment_count": fragment_count,
        "element_set": ",".join(elements),
        "has_salt": has_salt,
        "has_metal": has_metal,
        "has_isotope": has_isotope,
        "stereo_undefined": stereo_undefined,
        "pains_alert": (pains_count > 0) if pains_count is not None else "",
        "pains_alert_count": pains_count if pains_count is not None else "",
        "structure_flags": ";".join(flags),
        "morgan_radius": MORGAN_RADIUS,
        "morgan_nbits": MORGAN_NBITS,
        "morgan_use_chirality": MORGAN_USE_CHIRALITY,
        "morgan_fingerprint": fingerprint,
    }
    return values


def _index_external_rows(
    rows: Iterable[Mapping[str, object]],
    tool: str,
    *,
    master_by_id: Mapping[str, Mapping[str, object]],
) -> dict[str, dict[str, str]]:
    if tool not in TOOL_NAMES:
        raise CandidateValidationError(f"Unsupported external tool: {tool!r}")
    materialized = [dict(row) for row in rows]
    if not materialized:
        raise CandidateValidationError(f"{tool} result CSV is empty")
    indexed: dict[str, dict[str, str]] = {}
    for row_number, raw in enumerate(materialized, 2):
        row = {str(key): _text(value) for key, value in raw.items()}
        compound_id = row.get("compound_id", "")
        if not compound_id:
            raise CandidateValidationError(f"{tool} result row {row_number} has no compound_id")
        if compound_id in indexed:
            raise CandidateValidationError(f"Duplicate {tool} result compound_id: {compound_id!r}")
        if compound_id not in master_by_id:
            raise CandidateValidationError(f"Unknown {tool} result compound_id: {compound_id!r}")
        external_smiles = row.get("canonical_smiles") or row.get("smiles")
        expected_smiles = _text(master_by_id[compound_id].get("canonical_smiles"))
        if external_smiles and _text(external_smiles) != expected_smiles:
            raise CandidateValidationError(
                f"{tool} canonical_smiles mismatch for {compound_id!r}"
            )
        row["compound_id"] = compound_id
        row["canonical_smiles"] = expected_smiles
        indexed[compound_id] = row
    return indexed


def _external_result_state(row: Mapping[str, object] | None) -> str:
    if not row:
        return "missing"
    status = ""
    for key in ("result_status", "status", "prediction_status", "run_status"):
        if _text(row.get(key)):
            status = _text(row.get(key)).lower()
            break
    if status in _FAILURE_VALUES:
        return "failed"
    if status not in {"completed", "success", "present"}:
        return "pending"
    payload = row.get("prediction_payload_json", "")
    try:
        values = json.loads(str(payload))
    except (ValueError, TypeError):
        return "empty"
    if not isinstance(values, dict) or not values or not all(_text(v) for v in values.values()):
        return "empty"
    if any(str(v).strip().lower() in {"n/a", "na", "none", "null", "nan", "error", "failed", "pending", "-", "--"} for v in values.values()):
        return "failed"
    if not row.get("raw_file_sha256") or not row.get("raw_row_number"):
        return "unverified"
    return "present"



def _review_index(
    rows: Iterable[Mapping[str, object]] | None,
    *,
    master_by_id: Mapping[str, Mapping[str, object]],
) -> dict[str, dict[str, str]]:
    if rows is None:
        return {}
    materialized = [dict(row) for row in rows]
    indexed: dict[str, dict[str, str]] = {}
    for row_number, raw in enumerate(materialized, 2):
        row = {str(key): _text(value) for key, value in raw.items()}
        compound_id = row.get("compound_id", "")
        if not compound_id:
            raise CandidateValidationError(f"Review row {row_number} has no compound_id")
        if compound_id in indexed:
            raise CandidateValidationError(f"Duplicate human review compound_id: {compound_id!r}")
        if compound_id not in master_by_id:
            raise CandidateValidationError(f"Unknown human review compound_id: {compound_id!r}")
        if row.get("canonical_smiles") != master_by_id[compound_id]["canonical_smiles"]:
            raise CandidateValidationError("Human review structure mismatch: " + compound_id)
        indexed[compound_id] = row
    return indexed


def _review_decision(row: Mapping[str, object] | None) -> str:
    if not row:
        return "pending"
    decisions = set()
    for key in ("approval_status", "review_status", "decision", "approved"):
        value = _text(row.get(key)).lower()
        if not value: continue
        if value in {"approved", "approve", "yes", "true", "1", "pass", "passed"}:
            decisions.add("approved")
        elif value in {"rejected", "reject", "no", "false", "0", "fail", "failed"}:
            decisions.add("rejected")
        else:
            decisions.add("pending")
    if len(decisions)>1:
        raise CandidateValidationError("Conflicting human review fields")
    return next(iter(decisions), "pending")



def _qualification(
    *,
    descriptor: Mapping[str, object],
    tool_states: Mapping[str, str],
    review_decision: str,
) -> tuple[str, str]:
    blockers: list[str] = []
    if descriptor.get("rdkit_status") != "parsed":
        blockers.append(str(descriptor.get("rdkit_status") or "rdkit_unknown"))
    if not descriptor.get("morgan_fingerprint"):
        blockers.append("morgan_fingerprint_missing")
    for tool in TOOL_NAMES:
        state = tool_states.get(tool, "missing")
        if state != "present":
            blockers.append(f"{tool}_{state}")
    if review_decision == "rejected":
        return "rejected", "human_review_rejected"
    if review_decision != "approved":
        blockers.append("human_review_pending")
    if blockers:
        return "pending", ";".join(blockers)
    return "eligible", "three_tool_results_and_human_approval"


def assess_candidate_rows(
    master_rows: Iterable[Mapping[str, object]],
    ranking_rows_by_direction: Mapping[str, Iterable[Mapping[str, object]]] | None = None,
    *,
    external_rows_by_tool: Mapping[str, Iterable[Mapping[str, object]]] | None = None,
    review_rows: Iterable[Mapping[str, object]] | None = None,
) -> list[dict[str, object]]:
    """Build the lossless assessed master table.

    ``ranking_rows_by_direction`` may be omitted while assessment exports are
    being prepared.  If supplied, all ranking rows are validated and attached
    as ``rank_<direction>`` and ``drugclip_score_<direction>`` fields.
    """

    master = validate_master_rows(master_rows)
    master_by_id = {row["compound_id"]: row for row in master}
    rankings = validate_rankings(
        ranking_rows_by_direction or {}, master_ids=set(master_by_id)
    )
    external_rows_by_tool = external_rows_by_tool or {}
    external_indexes: dict[str, dict[str, dict[str, str]]] = {}
    for tool in TOOL_NAMES:
        raw_rows = external_rows_by_tool.get(tool)
        external_indexes[tool] = (
            _index_external_rows(raw_rows, tool, master_by_id=master_by_id)
            if raw_rows is not None
            else {}
        )
    review = _review_index(review_rows, master_by_id=master_by_id)

    assessed: list[dict[str, object]] = []
    for master_row in master:
        compound_id = master_row["compound_id"]
        output: dict[str, object] = dict(master_row)
        descriptor = assess_smiles(master_row["canonical_smiles"])
        output.update(descriptor)
        for direction in DIRECTIONS:
            output[f"rank_{direction}"] = ""
            output[f"drugclip_score_{direction}"] = ""
            for ranked in rankings.get(direction, []):
                if ranked["compound_id"] == compound_id:
                    output[f"rank_{direction}"] = ranked["rank"]
                    output[f"drugclip_score_{direction}"] = ranked["drugclip_score"]
                    break

        tool_states: dict[str, str] = {}
        for tool in TOOL_NAMES:
            tool_row = external_indexes[tool].get(compound_id)
            state = _external_result_state(tool_row)
            tool_states[tool] = state
            output[f"{tool}_result_state"] = state
            output[f"{tool}_raw_imported"] = state in {"present", "failed", "empty"}
            if tool_row:
                for key, value in tool_row.items():
                    if key in {"compound_id", "canonical_smiles", "smiles"}:
                        continue
                    output[f"{tool}__{key}"] = value
        review_row = review.get(compound_id)
        decision = _review_decision(review_row)
        evidence = {"compound_id": compound_id, "descriptor": descriptor,
                    "tools": {tool: external_indexes[tool].get(compound_id) for tool in TOOL_NAMES}}
        evidence_hash = hashlib.sha256(json.dumps(evidence, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()
        output["assessment_evidence_sha256"] = evidence_hash
        if decision == "approved":
            required_review = ("reviewer", "reviewed_at_utc", "assessment_comment", "assessment_evidence_sha256")
            if any(not _text((review_row or {}).get(k)) for k in required_review):
                raise CandidateValidationError("Approval is missing reviewer, date, comment or evidence: " + compound_id)
            if review_row["assessment_evidence_sha256"] != evidence_hash:
                raise CandidateValidationError("Review refers to changed assessment evidence: " + compound_id)
            if descriptor.get("structure_flags") and not _text(review_row.get("structure_flag_review")):
                raise CandidateValidationError("Flagged structure needs a review note: " + compound_id)
        output["human_approval_status"] = decision
        output["reviewer"] = _text(
            (review_row or {}).get("reviewer") or (review_row or {}).get("reviewer_name")
        )
        output["reviewed_at_utc"] = _text((review_row or {}).get("reviewed_at_utc"))
        output["assessment_comment"] = _text(
            (review_row or {}).get("assessment_comment")
            or (review_row or {}).get("review_note")
            or (review_row or {}).get("comment")
        )
        output["route_comment"] = _text((review_row or {}).get("route_comment"))
        state, reason = _qualification(
            descriptor=descriptor, tool_states=tool_states, review_decision=decision
        )
        output["assessment_state"] = state
        output["assessment_reason"] = reason
        output["eligible_for_selection"] = state == "eligible"
        output["scientific_status"] = "in_silico_assessment_only"
        assessed.append(output)
    return assessed


def tool_input_rows(master_rows: Iterable[Mapping[str, object]]) -> list[dict[str, str]]:
    master = validate_master_rows(master_rows)
    return [
        {"compound_id": row["compound_id"], "canonical_smiles": row["canonical_smiles"]}
        for row in master
    ]


def write_external_input_bundle(
    master_rows: Iterable[Mapping[str, object]], output_dir: Path
) -> dict[str, object]:
    """Write ordered SwissADME batches and full ADMETlab/ProTox inputs.

    SwissADME is split in input order into batches of at most 200 rows.  The
    manifest is itself evidence of order and batch boundaries; no external
    service is called by this function.
    """

    master = validate_master_rows(master_rows)
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = tool_input_rows(master)
    manifest: list[dict[str, object]] = []
    swiss_paths: list[str] = []
    for start in range(0, len(rows), MAX_SWISSADME_BATCH):
        batch_number = start // MAX_SWISSADME_BATCH + 1
        batch = rows[start : start + MAX_SWISSADME_BATCH]
        path = output_dir / f"swissadme_batch_{batch_number:03d}.csv"
        write_csv_rows(path, batch, ["compound_id", "canonical_smiles"])
        paste_path = output_dir / f"swissadme_batch_{batch_number:03d}.smi"
        paste_path.write_text("".join(r["canonical_smiles"]+" "+r["compound_id"]+"\n" for r in batch), encoding="utf-8")
        swiss_paths.append(str(paste_path))
        manifest.append(
            {
                "tool": "swissadme",
                "batch_number": batch_number,
                "path": str(path),
                "order_start": start + 1,
                "order_end": start + len(batch),
                "compound_count": len(batch),
                "sha256": sha256_file(path),
            }
        )
    for tool in ("admetlab3", "protox3"):
        path = output_dir / f"{tool}_input.csv"
        write_csv_rows(path, rows, ["compound_id", "canonical_smiles"])
        manifest.append(
            {
                "tool": tool,
                "batch_number": 1,
                "path": str(path),
                "order_start": 1,
                "order_end": len(rows),
                "compound_count": len(rows),
                "sha256": sha256_file(path),
            }
        )
    manifest_path = output_dir / "external_input_manifest.csv"
    write_csv_rows(
        manifest_path,
        manifest,
        ["tool", "batch_number", "path", "order_start", "order_end", "compound_count", "sha256"],
    )
    for tool in ("admetlab3", "protox3"):
        (output_dir / f"{tool}_input.smi").write_text("".join(r["canonical_smiles"]+" "+r["compound_id"]+"\n" for r in rows), encoding="utf-8")
    return {
        "compound_count": len(rows),
        "swissadme_batches": swiss_paths,
        "admetlab3_input": str(output_dir / "admetlab3_input.csv"),
        "protox3_input": str(output_dir / "protox3_input.csv"),
        "manifest": str(manifest_path),
        "manifest_sha256": sha256_file(manifest_path),
    }


REVIEW_TEMPLATE_FIELDS = (
    "compound_id",
    "canonical_smiles",
    "approval_status",
    "reviewer",
    "reviewed_at_utc",
    "assessment_comment",
    "route_comment",
    "structure_flag_review",
    "assessment_evidence_sha256",
)


def review_template_rows(master_rows: Iterable[Mapping[str, object]]) -> list[dict[str, object]]:
    master = validate_master_rows(master_rows)
    return [
        {
            "compound_id": row["compound_id"],
            "canonical_smiles": row["canonical_smiles"],
            "approval_status": "pending",
            "reviewer": "",
            "reviewed_at_utc": "",
            "assessment_comment": "",
            "route_comment": "",
            "structure_flag_review": "",
            "assessment_evidence_sha256": row.get("assessment_evidence_sha256", ""),
        }
        for row in master
    ]


def _fingerprint_bits(value: object, *, expected_length: int = MORGAN_NBITS) -> int:
    text = _text(value)
    if len(text) != expected_length or any(char not in "01" for char in text):
        raise CandidateValidationError(
            f"Fingerprint must be a {expected_length}-bit string for deterministic selection"
        )
    return int(text, 2)


def tanimoto_similarity(first: object, second: object) -> float:
    """Calculate exact Tanimoto similarity for two Morgan bit strings."""

    first_bits = _fingerprint_bits(first)
    second_bits = _fingerprint_bits(second)
    intersection = (first_bits & second_bits).bit_count()
    union = (first_bits | second_bits).bit_count()
    return 1.0 if union == 0 else intersection / union


def maxmin_select(
    rows: Iterable[Mapping[str, object]],
    count: int,
    *,
    excluded_ids: set[str] | None = None,
) -> list[dict[str, object]]:
    """Select a deterministic MaxMin subset from already eligible rows.

    The first seed is the best local rank.  Subsequent rows maximize their
    minimum Tanimoto distance to selected rows.  Exact ties use rank then
    compound_id, so input CSV order cannot change the result.
    """

    if count < 0:
        raise CandidateValidationError("Selection count must be non-negative")
    candidates: list[dict[str, object]] = []
    seen: set[str] = set()
    for row_number, raw in enumerate(rows, 1):
        row = dict(raw)
        compound_id = _text(row.get("compound_id"))
        if not compound_id:
            raise CandidateValidationError(f"Selection row {row_number} has no compound_id")
        if compound_id in seen:
            raise CandidateValidationError(f"Duplicate compound_id in selection pool: {compound_id!r}")
        rank = _safe_rank(row.get("rank"), row_number=row_number)
        _safe_float(row.get("drugclip_score", 0.0), field="drugclip_score", row_number=row_number)
        if not _text(row.get("morgan_fingerprint")):
            raise CandidateValidationError(f"Missing Morgan fingerprint for {compound_id!r}")
        _fingerprint_bits(row["morgan_fingerprint"])
        seen.add(compound_id)
        row["compound_id"] = compound_id
        row["rank"] = rank
        candidates.append(row)
    excluded_ids = excluded_ids or set()
    candidates = [row for row in candidates if str(row["compound_id"]) not in excluded_ids]
    candidates.sort(key=lambda row: (int(row["rank"]), str(row["compound_id"])))
    selected: list[dict[str, object]] = []
    remaining = list(candidates)
    while remaining and len(selected) < count:
        if not selected:
            chosen = remaining[0]
        else:
            scored: list[tuple[float, int, str, dict[str, object]]] = []
            for candidate in remaining:
                min_distance = min(
                    1.0
                    - tanimoto_similarity(
                        candidate["morgan_fingerprint"], selected_row["morgan_fingerprint"]
                    )
                    for selected_row in selected
                )
                scored.append(
                    (
                        min_distance,
                        int(candidate["rank"]),
                        str(candidate["compound_id"]),
                        candidate,
                    )
                )
            # Maximize distance, then prefer the local rank and stable ID.
            scored.sort(key=lambda item: (-item[0], item[1], item[2]))
            chosen = scored[0][3]
        selected.append(chosen)
        remaining.remove(chosen)
    return selected


def _ranked_selection_pool(
    assessed_rows: Iterable[Mapping[str, object]],
    direction: str,
    *,
    pool_size: int,
) -> list[dict[str, object]]:
    direction = _direction_from_row({}, direction) or ""
    if direction not in DIRECTIONS:
        raise CandidateValidationError(f"Unsupported direction: {direction!r}")
    pool: list[dict[str, object]] = []
    for row_number, raw in enumerate(assessed_rows, 2):
        row = dict(raw)
        if _normal_bool(row.get("eligible_for_selection")) != "true" or row.get("assessment_state") != "eligible":
            continue
        if row.get("human_approval_status") != "approved" or any(row.get(t+"_result_state") != "present" for t in TOOL_NAMES):
            raise CandidateValidationError("Inconsistent selection eligibility: " + str(row.get("compound_id")))
        rank_value = row.get(f"rank_{direction}")
        score_value = row.get(f"drugclip_score_{direction}")
        if not _text(rank_value):
            continue
        rank = _safe_rank(rank_value, row_number=row_number)
        score = _safe_float(score_value, field=f"drugclip_score_{direction}", row_number=row_number)
        row["rank"] = rank
        row["drugclip_score"] = score
        pool.append(row)
    pool.sort(key=lambda row: (int(row["rank"]), str(row["compound_id"])))
    # A ranking error inside the assessed table must stop selection instead
    # of being hidden by a top-N slice.
    ids = [str(row.get("compound_id", "")) for row in pool]
    if len(ids) != len(set(ids)):
        raise CandidateValidationError(f"Duplicate compound_id in eligible {direction} pool")
    ranks = [int(row["rank"]) for row in pool]
    if len(ranks) != len(set(ranks)):
        raise CandidateValidationError(f"Duplicate rank in eligible {direction} pool")
    return pool[:pool_size]


def select_candidates(
    assessed_rows: Iterable[Mapping[str, object]],
    *,
    directions: Sequence[str] = DIRECTIONS,
    pool_size: int = MAX_SELECTION_POOL,
    per_direction: int = DEFAULT_SELECTION_COUNT,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Select up to ``per_direction`` unique candidates per direction.

    Directions are processed in the declared order (the default is
    C006_C007, then C016).  Already selected IDs are excluded from later
    directions, allowing overlap to be resolved deterministically while
    preserving as many requested slots as the approved top-50 pools allow.
    """

    if pool_size < 1 or pool_size > MAX_SELECTION_POOL:
        raise CandidateValidationError(f"pool_size must be in 1..{MAX_SELECTION_POOL}")
    if per_direction < 1:
        raise CandidateValidationError("per_direction must be positive")
    normalized_directions: list[str] = []
    for direction in directions:
        normalized = _direction_from_row({}, direction) or ""
        if normalized not in DIRECTIONS:
            raise CandidateValidationError(f"Unsupported direction: {direction!r}")
        if normalized in normalized_directions:
            raise CandidateValidationError(f"Duplicate direction: {normalized!r}")
        normalized_directions.append(normalized)
    rows = [dict(row) for row in assessed_rows]
    for direction in normalized_directions:
        if rows and any(_text(r.get(f"rank_{direction}")) for r in rows):
            validate_ranking_rows([dict(compound_id=r.get("compound_id"),rank=r.get(f"rank_{direction}"),drugclip_score=r.get(f"drugclip_score_{direction}")) for r in rows],direction,master_ids={r["compound_id"] for r in rows})
    selected_ids: set[str] = set()
    selected_rows: list[dict[str, object]] = []
    summary: dict[str, object] = {
        "status": "complete",
        "requested_count": len(normalized_directions) * per_direction,
        "selected_count": 0,
        "directions": {},
        "parameters": {
            "directions": normalized_directions,
            "pool_size": pool_size,
            "per_direction": per_direction,
            "algorithm": "rank_seed_then_maxmin_tanimoto_distance",
            "morgan_radius": MORGAN_RADIUS,
            "morgan_nbits": MORGAN_NBITS,
            "morgan_use_chirality": MORGAN_USE_CHIRALITY,
            "tie_break": "local_rank_ascending_then_compound_id_ascending",
            "cross_direction_deduplication": "declared_direction_order",
        },
    }
    for direction in normalized_directions:
        pool = _ranked_selection_pool(rows, direction, pool_size=pool_size)
        available = [row for row in pool if str(row["compound_id"]) not in selected_ids]
        chosen = maxmin_select(available, per_direction)
        for position, row in enumerate(chosen, 1):
            item = dict(row)
            item["selection_direction"] = direction
            item["selection_position_in_direction"] = position
            item["selection_reason"] = (
                "best local rank seed"
                if position == 1
                else "max-min Morgan Tanimoto distance; ties resolved by rank then compound_id"
            )
            item["selection_parameters"] = json.dumps(
                summary["parameters"], ensure_ascii=False, sort_keys=True, separators=(",", ":")
            )
            item["source_rank"] = item["rank"]
            item["source_drugclip_score"] = item["drugclip_score"]
            selected_rows.append(item)
            selected_ids.add(str(item["compound_id"]))
        direction_summary = {
            "approved_ranked_pool_count": len(pool),
            "available_after_cross_direction_deduplication": len(available),
            "requested_count": per_direction,
            "selected_count": len(chosen),
            "status": "complete" if len(chosen) == per_direction else "partial",
            "overlap_excluded_count": len(pool) - len(available),
        }
        summary["directions"][direction] = direction_summary
        if len(chosen) < per_direction:
            summary["status"] = "partial"
    summary["selected_count"] = len(selected_rows)
    if len(selected_rows) < int(summary["requested_count"]):
        summary["status"] = "partial"
        summary["partial_reason"] = "approved non-duplicate candidates are insufficient"
    return selected_rows, summary


def canonical_assessment_fields(rows: Iterable[Mapping[str, object]]) -> list[str]:
    """Return stable columns with direction and tool provenance first."""

    materialized = [dict(row) for row in rows]
    preferred = [
        "compound_id",
        "canonical_smiles",
        "source",
        "source_ids",
        "source_tasks",
        "source_record_count",
        "rdkit_status",
        "rdkit_canonical_smiles",
        "molecular_weight",
        "mol_wt",
        "logp",
        "hbd",
        "hba",
        "rotatable_bonds",
        "tpsa",
        "heavy_atoms",
        "formal_charge",
        "ring_count",
        "fragment_count",
        "element_set",
        "has_salt",
        "has_metal",
        "has_isotope",
        "stereo_undefined",
        "pains_alert",
        "pains_alert_count",
        "structure_flags",
        "morgan_radius",
        "morgan_nbits",
        "morgan_use_chirality",
        "morgan_fingerprint",
    ]
    for direction in DIRECTIONS:
        preferred.extend([f"rank_{direction}", f"drugclip_score_{direction}"])
    for tool in TOOL_NAMES:
        preferred.append(f"{tool}_result_state")
    preferred.extend(
        [
            "human_approval_status",
            "reviewer",
            "reviewed_at_utc",
            "assessment_comment",
            "route_comment",
            "assessment_state",
            "assessment_reason",
            "eligible_for_selection",
            "scientific_status",
        ]
    )
    discovered: list[str] = []
    for field in preferred:
        if any(field in row for row in materialized) and field not in discovered:
            discovered.append(field)
    for row in materialized:
        for field in row:
            if field not in discovered:
                discovered.append(field)
    return discovered
