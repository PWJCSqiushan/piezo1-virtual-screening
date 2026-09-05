from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.pubchem import parse_property_csv


class PubChemTests(unittest.TestCase):
    def test_parse_property_csv_keeps_cid_and_smiles(self) -> None:
        rows = parse_property_csv(
            '"CID","SMILES","Title","MolecularFormula","MolecularWeight"\n'
            '2244,"CC(=O)OC1=CC=CC=C1C(=O)O","Aspirin","C9H8O4",180.16\n'
        )
        self.assertEqual(rows[0]["CID"], "2244")
        self.assertEqual(rows[0]["Title"], "Aspirin")

    def test_parse_property_csv_skips_missing_smiles(self) -> None:
        rows = parse_property_csv('"CID","SMILES"\n1,""\n2,"CCO"\n')
        self.assertEqual([row["CID"] for row in rows], ["2"])


if __name__ == "__main__":
    unittest.main()
