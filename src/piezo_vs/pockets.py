from __future__ import annotations

import csv
import math
import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

Residue = tuple[str, str]
Point3D = tuple[float, float, float]


@dataclass(frozen=True)
class Pocket:
    tool: str
    pocket_id: str
    center: Point3D
    residues: frozenset[Residue]
    source_file: str = ""
    rank: int | None = None
    score: float | None = None


def center_distance(left: Point3D, right: Point3D) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))


def _clean_row(row: dict[str, str]) -> dict[str, str]:
    return {(key or "").strip(): (value or "").strip() for key, value in row.items()}


def _parse_residue_token(token: str) -> Residue:
    chain, residue_number = token.split("_", 1)
    return chain or "_", residue_number


def parse_p2rank(path: Path, limit: int | None = None) -> list[Pocket]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [_clean_row(row) for row in csv.DictReader(handle)]
    pockets: list[Pocket] = []
    for row in rows[:limit]:
        residues = frozenset(_parse_residue_token(token) for token in row["residue_ids"].split())
        pockets.append(
            Pocket(
                tool="p2rank",
                pocket_id=row["name"],
                rank=int(row["rank"]),
                score=float(row["score"]),
                center=tuple(float(row[key]) for key in ("center_x", "center_y", "center_z")),
                residues=residues,
                source_file=str(path),
            )
        )
    return pockets


def parse_atom_pocket(path: Path, tool: str, pocket_id: str, rank: int | None = None) -> Pocket:
    residues: set[Residue] = set()
    coordinates: list[Point3D] = []
    with path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            chain = line[21].strip() or "_"
            residue_number = f"{line[22:26].strip()}{line[26].strip()}"
            try:
                point = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            residues.add((chain, residue_number))
            coordinates.append(point)
    if not coordinates:
        raise ValueError(f"No ATOM/HETATM coordinates found in {path}")
    count = len(coordinates)
    center = tuple(sum(point[index] for point in coordinates) / count for index in range(3))
    return Pocket(
        tool=tool,
        pocket_id=pocket_id,
        rank=rank,
        score=None,
        center=center,
        residues=frozenset(residues),
        source_file=str(path),
    )


def _natural_number(path: Path) -> int:
    numbers = re.findall(r"(\d+)", path.stem)
    return int(numbers[-1]) if numbers else 10**9


def parse_dogsite(download_dir: Path, pdb_id: str, limit: int | None = None) -> list[Pocket]:
    pattern = f"{pdb_id.lower()}_P_*_res.pdb"
    paths = sorted(download_dir.glob(pattern), key=_natural_number)
    if limit is not None:
        paths = paths[:limit]
    return [
        parse_atom_pocket(
            path,
            tool="dogsite3",
            pocket_id=path.name.removeprefix(f"{pdb_id.lower()}_").removesuffix("_res.pdb"),
            rank=index,
        )
        for index, path in enumerate(paths, 1)
    ]


def parse_fpocket(output_dir: Path, limit: int | None = None) -> list[Pocket]:
    pockets_dir = output_dir / "pockets" if (output_dir / "pockets").is_dir() else output_dir
    paths = sorted(pockets_dir.glob("pocket*_atm.pdb"), key=_natural_number)
    if limit is not None:
        paths = paths[:limit]
    if not paths:
        raise FileNotFoundError(f"No fpocket pocket*_atm.pdb files found under {pockets_dir}")
    return [
        parse_atom_pocket(
            path,
            tool="fpocket",
            pocket_id=path.stem.removesuffix("_atm"),
            rank=_natural_number(path),
        )
        for path in paths
    ]


def validate_fpocket_output_name(pdb_id: str, output_dir: Path) -> None:
    if output_dir.name.upper() != f"{pdb_id.upper()}_OUT":
        raise ValueError(
            f"fpocket output directory must be named {pdb_id.upper()}_out; got {output_dir.name}"
        )


def pockets_match(
    left: Pocket,
    right: Pocket,
    max_center_distance: float,
    min_shared_residues: int,
) -> bool:
    if left.tool == right.tool:
        return False
    return (
        center_distance(left.center, right.center) <= max_center_distance
        and len(left.residues & right.residues) >= min_shared_residues
    )


def build_consensus_components(
    pockets: Iterable[Pocket],
    max_center_distance: float = 12.0,
    min_shared_residues: int = 3,
) -> list[list[Pocket]]:
    items = list(pockets)
    by_tool: dict[str, list[Pocket]] = {}
    for pocket in items:
        by_tool.setdefault(pocket.tool, []).append(pocket)

    def identity(pocket: Pocket) -> tuple[str, str, str]:
        return pocket.tool, pocket.pocket_id, pocket.source_file

    def rank_sum(component: list[Pocket]) -> int:
        return sum(pocket.rank or 10**9 for pocket in component)

    def max_distance(component: list[Pocket]) -> float:
        return max(
            center_distance(left.center, right.center)
            for index, left in enumerate(component)
            for right in component[index + 1 :]
        )

    tool_names = sorted(by_tool)

    # A three-tool region is a strict clique: all three pairs must match and
    # at least min_shared_residues must be shared by all three tools. Source
    # pockets are allowed to participate in more than one hypothesis because
    # a broad pocket from one tool may legitimately overlap multiple regions
    # from another tool. Forcing one-to-one assignment would undercount tool
    # support and can falsely label an extendable pair as exact two-tool support.
    triples: list[list[Pocket]] = []
    if len(tool_names) >= 3:
        for first in by_tool[tool_names[0]]:
            for second in by_tool[tool_names[1]]:
                if not pockets_match(first, second, max_center_distance, min_shared_residues):
                    continue
                for third in by_tool[tool_names[2]]:
                    component = [first, second, third]
                    if not all(
                        pockets_match(left, right, max_center_distance, min_shared_residues)
                        for index, left in enumerate(component)
                        for right in component[index + 1 :]
                    ):
                        continue
                    if len(first.residues & second.residues & third.residues) < min_shared_residues:
                        continue
                    triples.append(component)
    triples.sort(
        key=lambda component: (
            -len(set.intersection(*(set(p.residues) for p in component))),
            max_distance(component),
            rank_sum(component),
        )
    )
    unique_triples: list[list[Pocket]] = []
    seen_components: set[frozenset[tuple[str, str, str]]] = set()
    for component in triples:
        identities = frozenset(identity(pocket) for pocket in component)
        if identities not in seen_components:
            unique_triples.append(component)
            seen_components.add(identities)
    triple_identity_sets = [
        {identity(pocket) for pocket in component} for component in unique_triples
    ]

    # A pair is exact two-tool support only when that same pair cannot be
    # extended by any third-tool pocket into a strict three-tool clique.
    pairs: list[list[Pocket]] = []
    for left_index, left_tool in enumerate(tool_names):
        for right_tool in tool_names[left_index + 1 :]:
            for left in by_tool[left_tool]:
                for right in by_tool[right_tool]:
                    if pockets_match(left, right, max_center_distance, min_shared_residues):
                        pairs.append([left, right])
    pairs.sort(
        key=lambda component: (
            -len(component[0].residues & component[1].residues),
            center_distance(component[0].center, component[1].center),
            rank_sum(component),
        )
    )
    exact_pairs: list[list[Pocket]] = []
    for component in pairs:
        identities = {identity(pocket) for pocket in component}
        if any(identities <= triple_identities for triple_identities in triple_identity_sets):
            continue
        frozen = frozenset(identities)
        if frozen not in seen_components:
            exact_pairs.append(component)
            seen_components.add(frozen)

    return sorted(
        [*unique_triples, *exact_pairs],
        key=lambda component: (
            -len(component),
            rank_sum(component),
            sorted((p.tool, p.pocket_id) for p in component),
        ),
    )


def group_overlapping_consensus_components(
    components: Iterable[list[Pocket]],
    *,
    max_center_distance: float = 12.0,
    min_shared_residues: int = 3,
) -> list[tuple[list[Pocket], list[list[Pocket]]]]:
    """Group redundant clique hypotheses around a stable best representative.

    Each candidate is a strict two- or three-tool clique. Candidates are sorted
    by support strength and residue overlap, then compared only with the fixed
    representative of each group. This avoids transitive chain merging while
    retaining an audit list of every alternative hypothesis.
    """

    def core(component: list[Pocket]) -> set[Residue]:
        return set.intersection(*(set(pocket.residues) for pocket in component))

    def center(component: list[Pocket]) -> Point3D:
        return tuple(
            sum(pocket.center[axis] for pocket in component) / len(component)
            for axis in range(3)
        )

    def spread(component: list[Pocket]) -> float:
        return max(
            (
                center_distance(left.center, right.center)
                for index, left in enumerate(component)
                for right in component[index + 1 :]
            ),
            default=0.0,
        )

    def rank_sum(component: list[Pocket]) -> int:
        return sum(pocket.rank or 10**9 for pocket in component)

    ordered = sorted(
        components,
        key=lambda component: (
            -len(component),
            -len(core(component)),
            spread(component),
            rank_sum(component),
            sorted((p.tool, p.pocket_id) for p in component),
        ),
    )
    groups: list[tuple[list[Pocket], list[list[Pocket]]]] = []
    for component in ordered:
        component_core = core(component)
        component_center = center(component)
        for representative, members in groups:
            if (
                center_distance(component_center, center(representative)) <= max_center_distance
                and len(component_core & core(representative)) >= min_shared_residues
            ):
                members.append(component)
                break
        else:
            groups.append((component, [component]))
    return groups


def summarize_component(component: list[Pocket]) -> dict[str, object]:
    tools = sorted({pocket.tool for pocket in component})
    tool_residues: dict[str, set[Residue]] = {tool: set() for tool in tools}
    for pocket in component:
        tool_residues[pocket.tool].update(pocket.residues)
    all_residues = set().union(*tool_residues.values())
    residue_support = {
        residue: sum(residue in residues for residues in tool_residues.values()) for residue in all_residues
    }
    core_two_plus = {residue for residue, support in residue_support.items() if support >= 2}
    core_all = {residue for residue, support in residue_support.items() if support == len(tools)}
    center = tuple(sum(p.center[index] for p in component) / len(component) for index in range(3))
    max_spread = max(
        (center_distance(left.center, right.center) for i, left in enumerate(component) for right in component[i + 1 :]),
        default=0.0,
    )
    pairwise_jaccard = []
    for left_index, left_tool in enumerate(tools):
        for right_tool in tools[left_index + 1 :]:
            left_residues = tool_residues[left_tool]
            right_residues = tool_residues[right_tool]
            union = left_residues | right_residues
            pairwise_jaccard.append(
                len(left_residues & right_residues) / len(union) if union else 0.0
            )
    minimum_jaccard_similarity = min(pairwise_jaccard, default=0.0)

    def residue_text(residues: set[Residue]) -> str:
        return " ".join(f"{chain}_{number}" for chain, number in sorted(residues))

    return {
        "support_count": len(tools),
        "tools": ",".join(tools),
        "member_count": len(component),
        "members": " ".join(f"{p.tool}:{p.pocket_id}" for p in sorted(component, key=lambda p: (p.tool, p.rank or 10**9))),
        "center_x": round(center[0], 3),
        "center_y": round(center[1], 3),
        "center_z": round(center[2], 3),
        "max_member_center_distance": round(max_spread, 3),
        "min_pairwise_jaccard_similarity": round(minimum_jaccard_similarity, 6),
        "max_pairwise_jaccard_distance": round(1.0 - minimum_jaccard_similarity, 6),
        "core_residue_count_2plus": len(core_two_plus),
        "core_residue_count_all_supporting_tools": len(core_all),
        "core_residues_2plus": residue_text(core_two_plus),
        "core_residues_all_supporting_tools": residue_text(core_all),
        "same_tool_duplicate_warning": any(sum(p.tool == tool for p in component) > 1 for tool in tools),
        "source_files": " | ".join(sorted({p.source_file for p in component})),
    }


def classify_pocket_scope(
    component: list[Pocket], piezo1_chains: set[str], mdfic_chains: set[str]
) -> tuple[str, str]:
    chains = {chain for pocket in component for chain, _ in pocket.residues}
    if chains and chains <= piezo1_chains:
        scope = "piezo1"
    elif mdfic_chains and chains <= mdfic_chains:
        scope = "mdfic"
    elif chains & piezo1_chains and chains & mdfic_chains:
        scope = "piezo1_mdfic_interface"
    else:
        scope = "other_or_unknown"
    return scope, ",".join(sorted(chains))
