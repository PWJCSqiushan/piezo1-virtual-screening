from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.drugclip_io import (
    extract_pocket_atoms,
    parse_ranked_compounds,
    parse_residue_tokens,
    read_compound_inputs,
    select_consensus_row,
)


class DrugClipIoTests(unittest.TestCase):
    def test_parse_residue_tokens(self) -> None:
        self.assertEqual(parse_residue_tokens("A_10 B_20A"), {("A", "10"), ("B", "20A")})

    def test_select_consensus_requires_latest_tier_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "consensus.csv"
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "pdb_id",
                        "consensus_id",
                        "classification_status",
                        "support_count",
                        "tier",
                    ],
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "pdb_id": "8YEZ",
                        "consensus_id": "C001",
                        "classification_status": "final_three_tool_run",
                        "support_count": "3",
                        "tier": "T2",
                    }
                )
            self.assertEqual(select_consensus_row(path, "8yez", "c001")["tier"], "T2")

    def test_extracts_only_selected_residue_and_primary_altloc(self) -> None:
        def atom_line(serial: int, atom: str, altloc: str, chain: str, residue: int, x: float) -> str:
            return (
                f"ATOM  {serial:5d} {atom:>4s}{altloc:1s}ALA {chain:1s}{residue:4d}    "
                f"{x:8.3f}{2.0:8.3f}{3.0:8.3f}{1.0:6.2f}{20.0:6.2f}          C \n"
            )

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.pdb"
            path.write_text(
                atom_line(1, "CA", "", "A", 10, 1.0)
                + atom_line(2, "CB", "B", "A", 10, 2.0)
                + atom_line(3, "CA", "", "B", 20, 3.0),
                encoding="ascii",
            )
            atoms, missing = extract_pocket_atoms(path, {("A", "10")})
            self.assertEqual([atom.atom_name for atom in atoms], ["CA"])
            self.assertFalse(missing)

    def test_reads_csv_compounds(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "molecules.csv"
            path.write_text("compound_id,smiles\naspirin,CC(=O)OC1=CC=CC=C1C(=O)O\n", encoding="utf-8")
            compounds = read_compound_inputs(path)
            self.assertEqual(compounds[0].compound_id, "aspirin")

    def test_parse_ranked_compounds(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ranked.txt"
            path.write_text("CCO\t0.75\nCCN\t-0.1\n", encoding="utf-8")
            self.assertEqual(parse_ranked_compounds(path), [("CCO", 0.75), ("CCN", -0.1)])


if __name__ == "__main__":
    unittest.main()

