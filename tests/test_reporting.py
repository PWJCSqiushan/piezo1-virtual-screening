from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.reporting import render_html


class ReportingTests(unittest.TestCase):
    def test_render_html_escapes_compound_fields(self) -> None:
        snapshot = {
            "active_pdb_ids": [],
            "total_consensus_regions": 0,
            "purpose": "technical_validation",
            "scientific_status": "computational_prediction_only",
            "generated_at_utc": "2026-09-06T00:00:00+00:00",
            "structures": [],
            "stages": [],
            "alerts": [],
            "drugclip": {
                "ranked_count": 1,
                "rows": [
                    {
                        "rank": "1",
                        "compound_id": "<unsafe>",
                        "source_id": "CID:1",
                        "source_url": "https://example.test/?x=1&y=2",
                        "drugclip_score": "0.5",
                    }
                ],
                "selected_consensus": {},
            },
            "docking": {"metadata": {}},
        }
        output = render_html(snapshot)
        self.assertIn("&lt;unsafe&gt;", output)
        self.assertNotIn(">\u003cunsafe>", output)
        self.assertIn("x=1&amp;y=2", output)


if __name__ == "__main__":
    unittest.main()
