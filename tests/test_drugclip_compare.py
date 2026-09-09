from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "19_compare_drugclip_runs.py"


def write_ranked(path: Path, rows: list[tuple[int, str, float]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["rank", "compound_id", "drugclip_score"]
        )
        writer.writeheader()
        for rank, compound_id, score in rows:
            writer.writerow(
                {"rank": rank, "compound_id": compound_id, "drugclip_score": score}
            )


class DrugClipComparisonTests(unittest.TestCase):
    def test_same_universe_is_comparable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.csv"
            second = root / "second.csv"
            output = root / "comparison.json"
            write_ranked(first, [(1, "A", 0.9), (2, "B", 0.8), (3, "C", 0.7)])
            write_ranked(second, [(1, "A", 0.91), (2, "B", 0.79), (3, "C", 0.69)])
            process = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--first",
                    str(first),
                    "--second",
                    str(second),
                    "--output",
                    str(output),
                    "--top-k",
                    "2",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(process.returncode, 0, process.stderr)
            result = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(result["comparison_status"], "comparable")
            self.assertEqual(result["spearman_rank_correlation"], 1.0)
            self.assertEqual(result["top_k_overlap_count"], 2)

    def test_different_universes_are_not_comparable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.csv"
            second = root / "second.csv"
            output = root / "comparison.json"
            write_ranked(first, [(1, "A", 0.9), (2, "B", 0.8)])
            write_ranked(second, [(1, "A", 0.9), (2, "C", 0.8)])
            process = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--first",
                    str(first),
                    "--second",
                    str(second),
                    "--output",
                    str(output),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(process.returncode, 2, process.stderr)
            result = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(result["comparison_status"], "not_comparable")
            self.assertIsNone(result["spearman_rank_correlation"])
            self.assertIn("compound_universes_differ", result["not_comparable_reasons"])


if __name__ == "__main__":
    unittest.main()