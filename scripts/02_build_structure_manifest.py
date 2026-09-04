from __future__ import annotations

import csv
import hashlib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.io_utils import nested_get, read_json, sha256_file, write_json
from piezo_vs.project import ensure_project_dirs

AA3 = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}


def join_values(value: object) -> str:
    if isinstance(value, list):
        return ",".join(str(item) for item in value)
    return "" if value is None else str(value)


def mutation_coordinate_status(pdb_path: Path, configured: dict[str, object]) -> tuple[str, int]:
    label = str(configured.get("mutation_label", "none"))
    if label.lower() in {"", "none", "wt"}:
        return "not_applicable_wild_type", 0
    substitution = re.fullmatch(r"([A-Z])(\d+)([A-Z])", label)
    deletion = re.fullmatch(r"([A-Z])(\d+)del", label, flags=re.IGNORECASE)
    match = substitution or deletion
    if not match or not pdb_path.exists():
        return "mutation_not_machine_parseable", 0
    position = int(match.group(2))
    piezo1_chains = set(configured.get("auth_piezo1_chains", []))
    residues: list[str] = []
    with pdb_path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  ") or line[21].strip() not in piezo1_chains:
                continue
            try:
                residue_number = int(line[22:26])
            except ValueError:
                continue
            if residue_number == position:
                residues.append(line[17:20].strip())
    if substitution:
        expected_residue = AA3[substitution.group(3)]
        if not residues:
            return "substitution_site_absent_from_coordinates", 0
        if expected_residue in residues:
            return f"substitution_modeled_as_{expected_residue}", len(residues)
        return f"substitution_site_present_but_not_{expected_residue}", len(residues)
    if residues:
        return "deletion_site_has_atoms_unexpected", len(residues)
    return "deletion_variant_not_directly_verifiable_from_atom_coordinates", 0


def coordinate_record_fingerprint(pdb_path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    count = 0
    with pdb_path.open("rb") as handle:
        for line in handle:
            if line.startswith((b"ATOM  ", b"HETATM")):
                digest.update(line.rstrip(b"\r\n"))
                digest.update(b"\n")
                count += 1
    return digest.hexdigest(), count


def main() -> int:
    ensure_project_dirs(ROOT)
    config = read_json(ROOT / "config" / "structures.json")
    rows: list[dict[str, object]] = []

    for configured in config["structures"]:
        pdb_id = configured["pdb_id"].upper()
        api_dir = ROOT / "data" / "raw" / "structures" / "api" / pdb_id
        entry_path = api_dir / "entry.json"
        if not entry_path.exists():
            raise FileNotFoundError(f"Missing {entry_path}; run 01_fetch_structures.py first")
        entry = read_json(entry_path)
        info = entry.get("rcsb_entry_info", {})

        entity_summaries: list[str] = []
        uniprots: set[str] = set()
        has_mdfic = False
        entity_ids = nested_get(entry, "rcsb_entry_container_identifiers", "polymer_entity_ids", default=[])
        for entity_id in entity_ids:
            entity_path = api_dir / f"polymer_entity_{entity_id}.json"
            if not entity_path.exists():
                continue
            entity = read_json(entity_path)
            description = nested_get(entity, "rcsb_polymer_entity", "pdbx_description")
            identifiers = entity.get("rcsb_polymer_entity_container_identifiers", {})
            auth_chains = join_values(identifiers.get("auth_asym_ids", []))
            asym_chains = join_values(identifiers.get("asym_ids", []))
            entity_uniprots = identifiers.get("uniprot_ids", []) or []
            uniprots.update(str(item) for item in entity_uniprots)
            description_upper = str(description).upper()
            has_mdfic = has_mdfic or "MDFIC" in description_upper or "Q9P1T7" in entity_uniprots
            entity_summaries.append(
                f"{entity_id}:{description};label_chains={asym_chains};auth_chains={auth_chains};"
                f"uniprot={join_values(entity_uniprots)}"
            )

        modeled = int(info.get("deposited_modeled_polymer_monomer_count", 0) or 0)
        total = int(info.get("deposited_polymer_monomer_count", 0) or 0)
        unmodeled = int(info.get("deposited_unmodeled_polymer_monomer_count", 0) or 0)
        cif_path = ROOT / "data" / "raw" / "structures" / f"{pdb_id}.cif"
        pdb_path = ROOT / "data" / "raw" / "structures" / f"{pdb_id}.pdb"
        coordinate_status, mutation_site_atom_count = mutation_coordinate_status(pdb_path, configured)
        coordinate_sha256, atom_record_count = coordinate_record_fingerprint(pdb_path)
        resolutions = info.get("resolution_combined", []) or []
        row = {
            "pdb_id": pdb_id,
            "role": configured["role"],
            "use_for_mvp": configured["use_for_mvp"],
            "category": configured.get("category", ""),
            "mutation_label": configured.get("mutation_label", ""),
            "mutation_coordinate_status": coordinate_status,
            "mutation_site_atom_count": mutation_site_atom_count,
            "coordinate_atom_sha256": coordinate_sha256,
            "atom_record_count": atom_record_count,
            "conformation_label": configured.get("conformation_label", ""),
            "expected_resolution_angstrom": configured.get("expected_resolution_angstrom", ""),
            "title": nested_get(entry, "struct", "title"),
            "experimental_method": info.get("experimental_method", ""),
            "resolution_angstrom": resolutions[0] if resolutions else "",
            "initial_release_date": nested_get(entry, "rcsb_accession_info", "initial_release_date"),
            "polymer_composition": info.get("polymer_composition", ""),
            "assembly_count": info.get("assembly_count", ""),
            "polymer_instances": info.get("deposited_polymer_entity_instance_count", ""),
            "modeled_residues": modeled,
            "unmodeled_residues": unmodeled,
            "deposited_residues": total,
            "modeled_fraction": round(modeled / total, 4) if total else "",
            "has_mdfic": has_mdfic,
            "uniprot_ids": ",".join(sorted(uniprots)),
            "polymer_entities": " | ".join(entity_summaries),
            "cif_sha256": sha256_file(cif_path) if cif_path.exists() else "",
            "pdb_sha256": sha256_file(pdb_path) if pdb_path.exists() else "",
            "notes": configured.get("notes", ""),
        }
        rows.append(row)

    csv_path = ROOT / "results" / "structure_manifest.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    write_json(ROOT / "results" / "structure_manifest.json", rows)
    print(f"Wrote {csv_path}")
    for row in rows:
        print(
            f"{row['pdb_id']}: resolution={row['resolution_angstrom']} angstrom, "
            f"modeled={row['modeled_fraction']}, MDFIC={row['has_mdfic']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
