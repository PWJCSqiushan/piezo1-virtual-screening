from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.pockets import (
    Pocket,
    build_consensus_components,
    center_distance,
    classify_pocket_scope,
    group_overlapping_consensus_components,
    summarize_component,
    validate_fpocket_output_name,
)


class ConsensusTests(unittest.TestCase):
    def pocket(self, tool: str, pocket_id: str, center: tuple[float, float, float], residues: set[tuple[str, str]]) -> Pocket:
        return Pocket(tool, pocket_id, center, frozenset(residues), source_file=f"{tool}.pdb", rank=1)

    def test_center_distance(self) -> None:
        self.assertEqual(center_distance((0, 0, 0), (3, 4, 0)), 5.0)

    def test_three_tool_region_has_support_three(self) -> None:
        common = {("A", "10"), ("A", "11"), ("A", "12")}
        pockets = [
            self.pocket("p2rank", "p1", (0, 0, 0), common | {("A", "13")}),
            self.pocket("dogsite3", "P_1", (1, 0, 0), common | {("A", "14")}),
            self.pocket("fpocket", "pocket1", (0, 1, 0), common | {("A", "15")}),
        ]
        components = build_consensus_components(pockets)
        self.assertEqual(len(components), 1)
        summary = summarize_component(components[0])
        self.assertEqual(summary["support_count"], 3)
        self.assertEqual(summary["core_residue_count_all_supporting_tools"], 3)

    def test_two_tool_region_is_kept(self) -> None:
        common = {("B", "20"), ("B", "21"), ("B", "22")}
        pockets = [
            self.pocket("p2rank", "p1", (0, 0, 0), common),
            self.pocket("dogsite3", "P_1", (2, 0, 0), common),
            self.pocket("fpocket", "pocket9", (50, 50, 50), common),
        ]
        components = build_consensus_components(pockets)
        self.assertEqual(len(components), 1)
        self.assertEqual(summarize_component(components[0])["support_count"], 2)

    def test_same_tool_predictions_do_not_form_consensus(self) -> None:
        common = {("A", "1"), ("A", "2"), ("A", "3")}
        pockets = [
            self.pocket("p2rank", "p1", (0, 0, 0), common),
            self.pocket("p2rank", "p2", (1, 0, 0), common),
        ]
        self.assertEqual(build_consensus_components(pockets), [])

    def test_transitive_chain_is_not_three_tool_consensus(self) -> None:
        ab = {("A", "1"), ("A", "2"), ("A", "3")}
        bc = {("A", "4"), ("A", "5"), ("A", "6")}
        pockets = [
            self.pocket("p2rank", "p1", (0, 0, 0), ab),
            self.pocket("dogsite3", "P_1", (10, 0, 0), ab | bc),
            self.pocket("fpocket", "pocket1", (20, 0, 0), bc),
        ]
        components = build_consensus_components(pockets)
        self.assertTrue(components)
        self.assertTrue(all(summarize_component(component)["support_count"] == 2 for component in components))

    def test_jaccard_similarity_and_distance_are_both_explicit(self) -> None:
        first = {("A", "1"), ("A", "2"), ("A", "3"), ("A", "4")}
        second = {("A", "1"), ("A", "2"), ("A", "3"), ("A", "5")}
        summary = summarize_component([
            self.pocket("p2rank", "p1", (0, 0, 0), first),
            self.pocket("fpocket", "f1", (1, 0, 0), second),
        ])
        self.assertAlmostEqual(summary["min_pairwise_jaccard_similarity"], 0.6)
        self.assertAlmostEqual(summary["max_pairwise_jaccard_distance"], 0.4)

    def test_three_pairs_without_triple_residue_overlap_is_not_support_three(self) -> None:
        ab = {("A", "1"), ("A", "2"), ("A", "3")}
        ac = {("A", "4"), ("A", "5"), ("A", "6")}
        bc = {("A", "7"), ("A", "8"), ("A", "9")}
        pockets = [
            self.pocket("p2rank", "p1", (0, 0, 0), ab | ac),
            self.pocket("dogsite3", "P_1", (1, 0, 0), ab | bc),
            self.pocket("fpocket", "pocket1", (0, 1, 0), ac | bc),
        ]
        components = build_consensus_components(pockets)
        self.assertTrue(all(summarize_component(component)["support_count"] == 2 for component in components))

    def test_pair_extendable_to_strict_triple_is_not_labeled_exact_two(self) -> None:
        common = {("A", "1"), ("A", "2"), ("A", "3")}
        pockets = [
            self.pocket("p2rank", "p1", (0, 0, 0), common),
            self.pocket("dogsite3", "d1", (1, 0, 0), common),
            self.pocket("fpocket", "f1", (0, 1, 0), common),
            self.pocket("dogsite3", "d2", (0, 2, 0), common | {("A", "4")}),
        ]
        components = build_consensus_components(pockets)
        extendable_pair = {("p2rank", "p1"), ("dogsite3", "d1")}
        pair_components = [
            {(pocket.tool, pocket.pocket_id) for pocket in component}
            for component in components
            if len(component) == 2
        ]
        self.assertNotIn(extendable_pair, pair_components)
        self.assertTrue(any(len(component) == 3 for component in components))

    def test_overlapping_pair_is_grouped_under_stronger_triple(self) -> None:
        common = {("A", "1"), ("A", "2"), ("A", "3")}
        triple = [
            self.pocket("p2rank", "p1", (0, 0, 0), common),
            self.pocket("dogsite3", "d1", (1, 0, 0), common),
            self.pocket("fpocket", "f1", (0, 1, 0), common),
        ]
        pair = [
            self.pocket("p2rank", "p2", (1, 1, 0), common | {("A", "4")}),
            self.pocket("dogsite3", "d2", (2, 1, 0), common | {("A", "5")}),
        ]
        groups = group_overlapping_consensus_components([pair, triple])
        self.assertEqual(len(groups), 1)
        representative, members = groups[0]
        self.assertEqual(len(representative), 3)
        self.assertEqual(len(members), 2)

    def test_fpocket_output_name_must_match_pdb(self) -> None:
        validate_fpocket_output_name("8YFC", Path("8YFC_out"))
        with self.assertRaises(ValueError):
            validate_fpocket_output_name("8YFC", Path("8ZU3_out"))

    def test_interface_scope_is_explicit(self) -> None:
        component = [
            self.pocket("p2rank", "p1", (0, 0, 0), {("A", "1"), ("C", "2")}),
            self.pocket("dogsite3", "P_1", (1, 0, 0), {("A", "1"), ("C", "2")}),
        ]
        scope, chains = classify_pocket_scope(component, {"A", "B", "D"}, {"C", "E", "F"})
        self.assertEqual(scope, "piezo1_mdfic_interface")
        self.assertEqual(chains, "A,C")


if __name__ == "__main__":
    unittest.main()
