from __future__ import annotations

import unittest
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.assessment import assessment_branch, index_tool_rows, prefixed_values


class AssessmentTests(unittest.TestCase):
    def test_oral_failure_is_retained_for_inhalation_review(self) -> None:
        branch, note = assessment_branch("No")
        self.assertEqual(branch, "inhalation_expert_review")
        self.assertIn("retained", note)

    def test_unknown_or_missing_value_needs_general_review(self) -> None:
        self.assertEqual(assessment_branch("")[0], "general_expert_review")
        self.assertEqual(assessment_branch("maybe")[0], "general_expert_review")

    def test_duplicate_tool_key_is_rejected(self) -> None:
        rows = [{"compound_id": "A"}, {"compound_id": "A"}]
        with self.assertRaises(ValueError):
            index_tool_rows(rows, "swissadme")

    def test_raw_columns_are_namespaced(self) -> None:
        values = prefixed_values({"compound_id": "A", "LogP": "2.1"}, "swissadme")
        self.assertEqual(values, {"swissadme__LogP": "2.1"})


if __name__ == "__main__":
    unittest.main()
