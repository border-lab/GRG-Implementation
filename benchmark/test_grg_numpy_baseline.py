"""Tests for grg_numpy_baseline (dense extraction + complementary NumPy recombination)."""

import numpy as np
import pytest

from grg_numpy_baseline import (
    estimate_numpy_memory,
    grg_to_numpy,
    grg_to_numpy_parallel,
    simulate_numpy_recombination,
)
from spot_check_numpy import complementary_pair_offspring


class FakeGRG:
    """Minimal GRG-like object (no PyGRGL)."""

    def __init__(self):
        self.num_nodes = 3
        self.num_samples = 2
        self.num_mutations = 2
        self._samples = {0, 1}
        self._up_edges = {0: [2], 1: [2], 2: []}
        self._node_mut_pairs = [(0, 0), (2, 1)]

    def is_sample(self, node_id):
        return node_id in self._samples

    def get_up_edges(self, node_id):
        return self._up_edges[node_id]

    def get_node_mutation_pairs(self):
        return list(self._node_mut_pairs)


class FakeBenchmarkComplementary:
    """One ``get_breakpoints`` result per mating pair (shared by sibling offspring)."""

    num_offspring_per_couple = 2

    def __init__(self, breakpoints_per_pair):
        self.breakpoints_per_pair = breakpoints_per_pair
        self.i = 0

    def get_breakpoints(self, bp_range, expected_crossovers=1.5):
        bp = np.asarray(self.breakpoints_per_pair[self.i], dtype=int)
        self.i += 1
        return bp, len(bp)


def test_estimate_numpy_memory_bytes_to_mb():
    mb = estimate_numpy_memory(100, 100, dtype_size=1)
    assert mb == pytest.approx(10000 / (1024 * 1024))


def test_grg_to_numpy_correct_matrix():
    g = FakeGRG()
    m = grg_to_numpy(g)
    assert m.shape == (2, 2)
    expected = np.array([[1, 1], [0, 1]], dtype=np.int8)
    np.testing.assert_array_equal(m, expected)


def test_grg_to_numpy_parallel_matches_serial():
    g = FakeGRG()
    serial = grg_to_numpy(g)
    parallel = grg_to_numpy_parallel(g, n_jobs=2)
    np.testing.assert_array_equal(parallel, serial)


def _assert_informative_columns_complement(parents, p1, p2, c0, c1):
    for col in range(parents.shape[1]):
        a1, a2 = parents[p1, col], parents[p2, col]
        v0, v1 = c0[col], c1[col]
        assert v0 in (a1, a2) and v1 in (a1, a2)
        if a1 != a2:
            assert v0 != v1


def test_simulate_numpy_recombination_complementary_one_pair():
    parents = np.array([[1, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 1]], dtype=np.int8)
    bench = FakeBenchmarkComplementary([[2, 4]])
    expected0, expected1 = complementary_pair_offspring(
        parents, 0, 1, np.array([2, 4]), first_child_starts_with_p1=True
    )

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(np.random, "shuffle", lambda x: None)
        mp.setattr(np.random, "randint", lambda low, high: 0)
        offspring, total_bp = simulate_numpy_recombination(
            bench, parents, bp_range=(0, 5), expected_crossovers=1.5
        )

    assert bench.i == 1
    assert total_bp == 2
    np.testing.assert_array_equal(offspring[0], expected0)
    np.testing.assert_array_equal(offspring[1], expected1)
    _assert_informative_columns_complement(parents, 0, 1, offspring[0], offspring[1])


def test_simulate_numpy_recombination_complementary_start_with_p2():
    parents = np.array([[1, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 1]], dtype=np.int8)
    bench = FakeBenchmarkComplementary([[2, 4]])
    expected0, expected1 = complementary_pair_offspring(
        parents, 0, 1, np.array([2, 4]), first_child_starts_with_p1=False
    )

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(np.random, "shuffle", lambda x: None)
        mp.setattr(np.random, "randint", lambda low, high: 1)
        offspring, total_bp = simulate_numpy_recombination(
            bench, parents, bp_range=(0, 5), expected_crossovers=1.5
        )

    assert total_bp == 2
    np.testing.assert_array_equal(offspring[0], expected0)
    np.testing.assert_array_equal(offspring[1], expected1)


def test_simulate_numpy_recombination_two_pairs():
    parents = np.array(
        [
            [1, 0, 1],
            [0, 1, 0],
            [1, 1, 0],
            [0, 0, 1],
        ],
        dtype=np.int8,
    )
    bench = FakeBenchmarkComplementary([[1], [2]])
    phase_seq = iter([1, 0])

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(np.random, "shuffle", lambda x: None)
        mp.setattr(np.random, "randint", lambda low, high: next(phase_seq))
        offspring, total_bp = simulate_numpy_recombination(
            bench, parents, bp_range=(0, 2), expected_crossovers=1.5
        )

    assert bench.i == 2
    assert total_bp == 2

    e00, e01 = complementary_pair_offspring(
        parents, 0, 1, np.array([1]), first_child_starts_with_p1=False
    )
    e10, e11 = complementary_pair_offspring(
        parents, 2, 3, np.array([2]), first_child_starts_with_p1=True
    )
    np.testing.assert_array_equal(offspring[0], e00)
    np.testing.assert_array_equal(offspring[1], e01)
    np.testing.assert_array_equal(offspring[2], e10)
    np.testing.assert_array_equal(offspring[3], e11)


def test_simulate_numpy_recombination_odd_num_samples():
    parents = np.array(
        [
            [1, 0, 1],
            [0, 1, 0],
            [1, 1, 1],
        ],
        dtype=np.int8,
    )
    bench = FakeBenchmarkComplementary([[]])

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(np.random, "shuffle", lambda x: None)
        mp.setattr(np.random, "randint", lambda low, high: 0)
        offspring, total_bp = simulate_numpy_recombination(
            bench, parents, bp_range=(0, 2), expected_crossovers=1.5
        )

    assert offspring.shape == (2, 3)
    assert total_bp == 0
    np.testing.assert_array_equal(offspring[0], parents[0])
    np.testing.assert_array_equal(offspring[1], parents[1])


class FakeBenchmarkK:
    """Shared breakpoints per pair; ``num_offspring_per_couple`` = k."""

    def __init__(self, k, breakpoints_per_pair):
        self.num_offspring_per_couple = k
        self.breakpoints_per_pair = breakpoints_per_pair
        self.i = 0

    def get_breakpoints(self, bp_range, expected_crossovers=1.5):
        bp = np.asarray(self.breakpoints_per_pair[self.i], dtype=int)
        self.i += 1
        return bp, len(bp)


def test_simulate_numpy_three_offspring_shared_breakpoints():
    """k==3: one breakpoint draw per pair; siblings partition p1/p2 per segment."""
    parents = np.array([[1, 0], [0, 1]], dtype=np.int8)
    bench = FakeBenchmarkK(3, [[1]])
    expected = np.array([[1, 1], [0, 1], [0, 0]], dtype=np.int8)

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(np.random, "shuffle", lambda x: None)
        mp.setattr(np.random, "randint", lambda low, high: 0)
        offspring, total_bp = simulate_numpy_recombination(
            bench, parents, bp_range=(0, 1), expected_crossovers=1.5
        )

    assert offspring.shape == (3, 2)
    assert total_bp == 1
    assert bench.i == 1
    np.testing.assert_array_equal(offspring, expected)
