from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.docking import docking_box_from_atoms, write_receptor_pdb
from piezo_vs.drugclip_io import PocketAtom


class DockingTests(unittest.TestCase):
    def test_box_has_minimum_size_and_reports_clipping(self) -> None:
        atoms = [
            PocketAtom("C", (-30.0, 0.0, 0.0), "A", "1"),
            PocketAtom("C", (30.0, 1.0, 1.0), "A", "2"),
        ]
        box = docking_box_from_atoms(atoms, (0.0, 0.0, 0.0))
        self.assertEqual(box["size_y"], 18.0)
        self.assertEqual(box["size_x"], 40.0)
        self.assertTrue(box["clipped_to_maximum"])

    def test_receptor_export_keeps_only_selected_chain_and_altloc(self) -> None:
        def line(serial: int, altloc: str, chain: str) -> str:
            return (
                f"ATOM  {serial:5d} {'CA':>4s}{altloc:1s}ALA {chain:1s}{10:4d}    "
                f"{1.0:8.3f}{2.0:8.3f}{3.0:8.3f}{1.0:6.2f}{20.0:6.2f}          C \n"
            )

        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "source.pdb"
            output = Path(directory) / "receptor.pdb"
            source.write_text(line(1, "", "A") + line(2, "B", "A") + line(3, "", "C"), encoding="ascii")
            count = write_receptor_pdb(source, output, {"A"})
            self.assertEqual(count, 1)
            self.assertIn("PRIMARY ALTLOC", output.read_text(encoding="ascii"))


if __name__ == "__main__":
    unittest.main()
