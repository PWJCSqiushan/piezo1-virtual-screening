from __future__ import annotations

import csv
import pickle
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


@dataclass(frozen=True)
class CompoundInput:
    compound_id: str
    smiles: str
    source_row: int


@dataclass(frozen=True)
class PocketAtom:
    atom_name: str
    coordinate: tuple[float, float, float]
    chain: str
    residue_number: str


def parse_residue_tokens(text: str) -> set[tuple[str, str]]:
    residues: set[tuple[str, str]] = set()
    for token in text.split():
        if "_" not in token:
            raise ValueError(f"Invalid consensus residue token: {token!r}")
        chain, residue_number = token.split("_", 1)
        if not chain or not residue_number:
            raise ValueError(f"Invalid consensus residue token: {token!r}")
        residues.add((chain, residue_number))
    if not residues:
        raise ValueError("The selected consensus pocket has no shared residues")
    return residues


def select_consensus_row(path: Path, pdb_id: str, consensus_id: str) -> dict[str, str]:
    pdb_id = pdb_id.upper()
    consensus_id = consensus_id.upper()
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    matches = [
        row
        for row in rows
        if (row.get("pdb_id") or "").upper() == pdb_id
        and (row.get("consensus_id") or "").upper() == consensus_id
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one row for {pdb_id}/{consensus_id} in {path}; found {len(matches)}"
        )
    row = matches[0]
    if row.get("classification_status") != "final_three_tool_run":
        raise ValueError("DrugCLIP input requires a final three-tool consensus run")
    support_count = int(row.get("support_count") or 0)
    if support_count not in {2, 3}:
        raise ValueError(f"Unsupported consensus support_count={support_count}")
    expected_tier = {2: "T1", 3: "T2"}[support_count]
    if row.get("tier") != expected_tier:
        raise ValueError(
            f"Tier mismatch for {pdb_id}/{consensus_id}: support={support_count} "
            f"requires {expected_tier}, got {row.get('tier')!r}"
        )
    return row


def extract_pocket_atoms(
    pdb_path: Path, residues: set[tuple[str, str]]
) -> tuple[list[PocketAtom], set[tuple[str, str]]]:
    atoms: list[PocketAtom] = []
    found_residues: set[tuple[str, str]] = set()
    with pdb_path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  "):
                continue
            altloc = line[16].strip()
            if altloc not in {"", "A"}:
                continue
            chain = line[21].strip() or "_"
            residue_number = f"{line[22:26].strip()}{line[26].strip()}"
            residue = (chain, residue_number)
            if residue not in residues:
                continue
            try:
                coordinate = (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
            except ValueError:
                continue
            atom_name = line[12:16].strip()
            if not atom_name:
                continue
            atoms.append(PocketAtom(atom_name, coordinate, chain, residue_number))
            found_residues.add(residue)
    return atoms, residues - found_residues


def _header_lookup(fieldnames: Iterable[str | None]) -> dict[str, str]:
    return {(name or "").strip().lower(): name or "" for name in fieldnames}


def read_compound_inputs(path: Path) -> list[CompoundInput]:
    suffix = path.suffix.lower()
    if suffix in {".csv", ".tsv"}:
        delimiter = "\t" if suffix == ".tsv" else ","
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            lookup = _header_lookup(reader.fieldnames or [])
            smiles_key = next(
                (lookup[key] for key in ("smiles", "smi", "canonical_smiles") if key in lookup),
                None,
            )
            id_key = next(
                (lookup[key] for key in ("compound_id", "molecule_id", "id", "name") if key in lookup),
                None,
            )
            if not smiles_key:
                raise ValueError(f"{path} needs a SMILES/smi/canonical_smiles column")
            compounds = []
            for row_number, row in enumerate(reader, 2):
                smiles = (row.get(smiles_key) or "").strip()
                if not smiles:
                    continue
                compound_id = (row.get(id_key) or "").strip() if id_key else ""
                compounds.append(CompoundInput(compound_id or f"row_{row_number}", smiles, row_number))
            return compounds

    compounds: list[CompoundInput] = []
    with path.open(encoding="utf-8-sig") as handle:
        for row_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            if row_number == 1 and parts[0].lower() in {"smiles", "smi", "canonical_smiles"}:
                continue
            compounds.append(
                CompoundInput(parts[1] if len(parts) > 1 else f"row_{row_number}", parts[0], row_number)
            )
    return compounds


def write_lmdb_records(path: Path, records: list[dict[str, Any]]) -> None:
    try:
        import lmdb
    except ImportError as exc:
        raise RuntimeError("Missing dependency 'lmdb'; run the DrugCLIP bootstrap first") from exc
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing LMDB: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    payloads = [pickle.dumps(record, protocol=pickle.HIGHEST_PROTOCOL) for record in records]
    estimated = sum(len(payload) for payload in payloads)
    map_size = max(64 * 1024 * 1024, estimated * 8)
    environment = lmdb.open(
        str(path),
        subdir=False,
        readonly=False,
        lock=False,
        readahead=False,
        meminit=False,
        map_size=map_size,
    )
    try:
        with environment.begin(write=True) as transaction:
            for index, payload in enumerate(payloads):
                transaction.put(str(index).encode("ascii"), payload)
        environment.sync()
    finally:
        environment.close()


def inspect_lmdb(path: Path) -> tuple[int, dict[str, Any]]:
    try:
        import lmdb
    except ImportError as exc:
        raise RuntimeError("Missing dependency 'lmdb'; run the DrugCLIP bootstrap first") from exc
    environment = lmdb.open(
        str(path), subdir=False, readonly=True, lock=False, readahead=False, meminit=False
    )
    try:
        with environment.begin() as transaction:
            count = int(transaction.stat()["entries"])
            cursor = transaction.cursor()
            first = next(iter(cursor), None)
            if first is None:
                return 0, {}
            return count, pickle.loads(first[1])
    finally:
        environment.close()


def parse_ranked_compounds(path: Path) -> list[tuple[str, float]]:
    rows: list[tuple[str, float]] = []
    with path.open(encoding="utf-8") as handle:
        for row_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            try:
                smiles, score_text = line.rsplit("\t", 1)
                score = float(score_text)
            except ValueError as exc:
                raise ValueError(f"Invalid DrugCLIP result at line {row_number}: {line!r}") from exc
            rows.append((smiles, score))
    return rows

