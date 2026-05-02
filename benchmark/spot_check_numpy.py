"""
Spot-check for *complementary* recombination (meiosis-style partitioning).

For each mating pair we use **one** crossover schedule (shared breakpoints). The two
offspring copy **opposite** parents in every segment: wherever child A takes a segment
from p1, child B takes that segment from p2, and vice versa. No genomic interval is
inherited from the same parent by both siblings at once.

This matches the intended biology for two products of one crossover process.
``simulate_numpy_recombination`` in grg_numpy_baseline.py uses the same rule when
``num_offspring_per_couple`` is 2.

Run from GRG-Implementation/benchmark::

    python spot_check_recombination.py
"""

from __future__ import annotations

import numpy as np


def complementary_pair_offspring(
    genotype_matrix: np.ndarray,
    p1: int,
    p2: int,
    bp: np.ndarray,
    first_child_starts_with_p1: bool,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build two offspring rows for one (p1, p2) pair with one breakpoint array ``bp``.

    Child 0 alternates parents along segments starting with p1 or p2 according to
    ``first_child_starts_with_p1``. Child 1 always uses the other parent in each
    corresponding segment.
    """
    genome_length = genotype_matrix.shape[1]
    bp = np.asarray(bp, dtype=int)
    out0 = np.empty(genome_length, dtype=genotype_matrix.dtype)
    out1 = np.empty(genome_length, dtype=genotype_matrix.dtype)

    if first_child_starts_with_p1:
        cur0, cur1 = p1, p2
    else:
        cur0, cur1 = p2, p1

    start_idx = 0
    for b in bp:
        out0[start_idx:b] = genotype_matrix[cur0, start_idx:b]
        out1[start_idx:b] = genotype_matrix[cur1, start_idx:b]
        start_idx = b
        cur0, cur1 = cur1, cur0

    out0[start_idx:] = genotype_matrix[cur0, start_idx:]
    out1[start_idx:] = genotype_matrix[cur1, start_idx:]
    return out0, out1


def _parent_label(row_idx: int, p1: int, p2: int) -> str:
    if row_idx == p1:
        return "p1"
    if row_idx == p2:
        return "p2"
    return str(row_idx)


def trace_complementary_pair(
    genotype_matrix: np.ndarray,
    p1: int,
    p2: int,
    bp: np.ndarray,
    first_child_starts_with_p1: bool,
) -> tuple[np.ndarray, np.ndarray, str]:
    """Same as complementary_pair_offspring but returns a human-readable trace."""
    genome_length = genotype_matrix.shape[1]
    bp = np.asarray(bp, dtype=int)
    lines: list[str] = []

    if first_child_starts_with_p1:
        cur0, cur1 = p1, p2
        lines.append(
            "  Child 0 starts on p1; child 1 starts on p2 (segment-wise complement)."
        )
    else:
        cur0, cur1 = p2, p1
        lines.append(
            "  Child 0 starts on p2; child 1 starts on p1 (segment-wise complement)."
        )

    out0 = np.empty(genome_length, dtype=genotype_matrix.dtype)
    out1 = np.empty(genome_length, dtype=genotype_matrix.dtype)

    start_idx = 0
    seg = 1
    for b in bp:
        out0[start_idx:b] = genotype_matrix[cur0, start_idx:b]
        out1[start_idx:b] = genotype_matrix[cur1, start_idx:b]
        lines.append(
            f"  Segment {seg}: cols [{start_idx}:{b}) — "
            f"child0 ← {_parent_label(cur0, p1, p2)} {np.asarray(out0[start_idx:b]).tolist()}, "
            f"child1 ← {_parent_label(cur1, p1, p2)} {np.asarray(out1[start_idx:b]).tolist()}"
        )
        start_idx = b
        cur0, cur1 = cur1, cur0
        seg += 1

    out0[start_idx:] = genotype_matrix[cur0, start_idx:]
    out1[start_idx:] = genotype_matrix[cur1, start_idx:]
    lines.append(
        f"  Final segment: cols [{start_idx}:{genome_length}] — "
        f"child0 ← {_parent_label(cur0, p1, p2)} {np.asarray(out0[start_idx:]).tolist()}, "
        f"child1 ← {_parent_label(cur1, p1, p2)} {np.asarray(out1[start_idx:]).tolist()}"
    )
    lines.append(f"  Child 0 full row: {out0.tolist()}")
    lines.append(f"  Child 1 full row: {out1.tolist()}")

    return out0, out1, "\n".join(lines)


def assert_source_partition_per_column(
    genotype_matrix: np.ndarray,
    p1: int,
    p2: int,
    child0: np.ndarray,
    child1: np.ndarray,
) -> None:
    """
    For every column where parental alleles differ, the two offspring alleles must
    differ (each copies a different parent at that locus under complementarity).
    """
    for col in range(genotype_matrix.shape[1]):
        a1, a2 = genotype_matrix[p1, col], genotype_matrix[p2, col]
        v0, v1 = child0[col], child1[col]
        assert v0 in (a1, a2) and v1 in (a1, a2)
        if a1 != a2:
            assert v0 != v1, (
                f"col {col}: parental alleles differ ({a1} vs {a2}) but both offspring "
                f"have allele {v0}; expected complementary haplotype sources."
            )


def spot_check_complementary(
    parents_2d: np.ndarray,
    breakpoints_per_pair: list[list[int]],
    first_child_starts_with_p1_per_pair: list[bool],
) -> np.ndarray:
    """
    Build offspring matrix using complementary siblings only.

    Parameters
    ----------
    parents_2d
        Rows are samples; consecutive pairs (0,1), (2,3), … mate.
    breakpoints_per_pair
        One list of breakpoint column indices per mating pair (not per offspring).
    first_child_starts_with_p1_per_pair
        One coin flip per pair: whether child 0's first segment comes from p1.
    """
    num_samples = parents_2d.shape[0]
    num_pairs = num_samples // 2
    assert len(breakpoints_per_pair) == num_pairs
    assert len(first_child_starts_with_p1_per_pair) == num_pairs

    offspring = np.zeros((num_pairs * 2, parents_2d.shape[1]), dtype=parents_2d.dtype)

    parent_indices = np.arange(num_samples)
    child_row = 0
    for pair_k in range(num_pairs):
        i = 2 * pair_k
        p1 = int(parent_indices[i])
        p2 = int(parent_indices[i + 1])
        bp = np.asarray(breakpoints_per_pair[pair_k], dtype=int)
        start_p1 = first_child_starts_with_p1_per_pair[pair_k]

        print(f"\n=== Pair {pair_k}: parent rows p1={p1}, p2={p2} ===")
        print(f"  Shared breakpoints (one crossover schedule): {bp.tolist()}")
        print(f"  Parent p1: {parents_2d[p1].tolist()}")
        print(f"  Parent p2: {parents_2d[p2].tolist()}")

        c0, c1, log = trace_complementary_pair(parents_2d, p1, p2, bp, start_p1)
        print(log)

        offspring[child_row] = c0
        offspring[child_row + 1] = c1

        assert_source_partition_per_column(parents_2d, p1, p2, c0, c1)

        child_row += 2

    print("\nComplementary spot check passed (distinct alleles at informative SNPs).")
    return offspring


def compare_to_current_simulate_numpy(
    parents_2d: np.ndarray,
    breakpoints_per_pair: list[list[int]],
    first_child_starts_with_p1_per_pair: list[bool],
    bp_range: tuple[int, int],
) -> None:
    """
    Run ``simulate_numpy_recombination`` with deterministic shuffle and phase
    (``numpy.random.randint``) for comparison to the manual trace above.
    """
    from unittest import mock

    from grg_numpy_baseline import simulate_numpy_recombination

    num_pairs = parents_2d.shape[0] // 2

    class Bench:
        num_offspring_per_couple = 2

        def __init__(self, per_pair):
            self.per_pair = per_pair
            self.i = 0

        def get_breakpoints(self, bp_range_inner, expected_crossovers=1.5):
            bp = np.asarray(self.per_pair[self.i], dtype=int)
            self.i += 1
            return bp, len(bp)

    # Phase for slot rotation (k==2): 0 matches first_child_starts_with_p1=True, 1 matches False.
    phases = [0 if start_p1 else 1 for start_p1 in first_child_starts_with_p1_per_pair]

    bench = Bench(breakpoints_per_pair)
    with mock.patch("numpy.random.shuffle", lambda x: None), mock.patch(
        "numpy.random.randint", side_effect=phases
    ):
        sim_offspring, total_bp = simulate_numpy_recombination(
            bench, parents_2d, bp_range=bp_range, expected_crossovers=1.5
        )

    expected_bp_count = sum(len(breakpoints_per_pair[k]) for k in range(num_pairs))
    print("\n--- simulate_numpy_recombination (matches complementary spot check) ---")
    print(f"  total_bp: {total_bp} (expected {expected_bp_count} for this schedule).")
    print("  Offspring matrix:\n", sim_offspring)


if __name__ == "__main__":
    parents = np.array(
        [
            [1, 1, 1, 0, 0, 0],
            [0, 0, 0, 1, 1, 1],
        ],
        dtype=np.int8,
    )

    offspring = spot_check_complementary(
        parents,
        breakpoints_per_pair=[[2, 4]],
        first_child_starts_with_p1_per_pair=[True],
    )
    print("\nFinal complementary offspring matrix (rows = offspring):")
    print(offspring)

    compare_to_current_simulate_numpy(
        parents,
        breakpoints_per_pair=[[2, 4]],
        first_child_starts_with_p1_per_pair=[True],
        bp_range=(0, 5),
    )
