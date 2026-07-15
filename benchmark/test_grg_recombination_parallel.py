"""Correctness tests for `ParallelNonDuplicationRecombination`.

Compares the parallel (node-aggregated, discover+commit) algorithm against:
  1. `NonDuplicationRecombination` (the sequential reference), on identical
     input graphs and identical recombination segments.
  2. A ground truth computed directly from graph reachability + the
     requested segments, independent of either recombination algorithm.

No PyGRGL dependency: `FakeGRG` is a minimal duck-typed stand-in exposing
exactly the surface `grg_recombination(_parallel).py` calls, following the
same pattern as `FakeGRG` in `test_grg_numpy_baseline.py`. `make_node`'s
`negative` semantics mirror upstream `grgl::MutableGRG::makeNode` (see
`grgl/include/grgl/grg.h`): the returned `NodeID` is *always* the plain,
non-negative slot index -- `negative=True` only affects internal
topological-order bookkeeping upstream (not modeled here). It is the
caller's responsibility to negate the id before passing it to `connect()`
whenever it should be treated as a "negative" node (this is exactly what
`NonDuplicationRecombination._recurse_attach`'s `neg_offspring = -offspring_id`
and `_extract_bubble`'s `connect(bubble_id, -offspring_id)` do). `_resolve`
mirrors upstream `nodeStripNegative`: sign-flip if negative, else identity.
"""

import random
from collections import deque

import pytest

from grg_recombination import NonDuplicationRecombination, simulate_grg_recombination
from grg_recombination_parallel import (
    ParallelNonDuplicationRecombination,
    _interval_covers_range,
    _interval_diff,
    _interval_overlaps_range,
    _interval_union,
    _normalize_intervals,
)


# ---------------------------------------------------------------------------
# Interval-set algebra unit tests (pure functions, no graph needed)
# ---------------------------------------------------------------------------

def _to_set(intervals):
    """Expand a canonical interval-set into a Python set of integer points,
    for brute-force comparison over a bounded universe."""
    s = set()
    for lo, hi in intervals:
        s.update(range(lo, hi))
    return s


def _random_intervals(rng, universe=40, n=4):
    pts = sorted(rng.sample(range(universe), 2 * n))
    return [(pts[i], pts[i + 1]) for i in range(0, len(pts), 2) if pts[i] < pts[i + 1]]


def test_normalize_merges_overlaps_and_adjacent():
    assert _normalize_intervals([(0, 3), (3, 5), (10, 12), (1, 2)]) == [(0, 5), (10, 12)]
    assert _normalize_intervals([(5, 5), (1, 1), (2, 4)]) == [(2, 4)]
    assert _normalize_intervals([]) == []


def test_interval_diff_matches_bruteforce():
    rng = random.Random(0)
    for _ in range(200):
        a = _normalize_intervals(_random_intervals(rng))
        b = _normalize_intervals(_random_intervals(rng))
        got = _interval_diff(a, b)
        assert got == _normalize_intervals(got), "diff output must already be canonical"
        assert _to_set(got) == _to_set(a) - _to_set(b)


def test_interval_union_matches_bruteforce():
    rng = random.Random(1)
    for _ in range(200):
        a = _normalize_intervals(_random_intervals(rng))
        b = _normalize_intervals(_random_intervals(rng))
        got = _interval_union(a, b)
        assert got == _normalize_intervals(got)
        assert _to_set(got) == _to_set(a) | _to_set(b)


def test_interval_overlaps_and_covers_range():
    rng = random.Random(2)
    for _ in range(200):
        ivs = _normalize_intervals(_random_intervals(rng))
        lo, hi = sorted(rng.sample(range(40), 2))
        if lo == hi:
            continue
        got_overlap = _interval_overlaps_range(ivs, lo, hi)
        expected_overlap = bool(_to_set(ivs) & set(range(lo, hi)))
        assert got_overlap == expected_overlap

        got_cover = _interval_covers_range(ivs, lo, hi)
        expected_cover = set(range(lo, hi)) <= _to_set(ivs)
        assert got_cover == expected_cover


# ---------------------------------------------------------------------------
# Minimal duck-typed GRG mock
# ---------------------------------------------------------------------------

class _Mutation:
    __slots__ = ("id", "position")

    def __init__(self, id_, position):
        self.id = id_
        self.position = position


class FakeGRG:
    """Minimal mutable-GRG-like object covering exactly the API surface used
    by `NonDuplicationRecombination` / `ParallelNonDuplicationRecombination`.

    Negative-id convention (matches real `pygrgl`/`grgl`): `make_node(negative=True)`
    returns the plain non-negative slot id it just allocated -- the `negative`
    flag does not affect the returned value's sign at all upstream. Callers
    that want a node to behave as "negative" (topologically-prepended) must
    negate the id themselves at every `connect()` call site (see
    `NonDuplicationRecombination._recurse_attach`'s `neg_offspring = -offspring_id`).
    `_resolve` mirrors `nodeStripNegative`: sign-flip if negative, else identity.
    """

    def __init__(self, bp_range):
        self.bp_range = bp_range
        self._up = {}
        self._down = {}
        self._muts_by_node = {}
        self._mutations = {}
        self._next_node_id = 0
        self._next_mut_id = 0
        self._samples = set()

    def _resolve(self, node_id):
        return -node_id if node_id < 0 else node_id

    @property
    def num_nodes(self):
        return self._next_node_id

    @property
    def num_edges(self):
        return sum(len(v) for v in self._down.values())

    @property
    def num_mutations(self):
        return len(self._mutations)

    @property
    def num_samples(self):
        return len(self._samples)

    def is_sample(self, node_id):
        return self._resolve(node_id) in self._samples

    def get_root_nodes(self):
        return [n for n in range(self._next_node_id) if not self._up.get(n)]

    def get_sample_nodes(self):
        return list(self._samples)

    def set_samples(self, ids):
        self._samples = set(self._resolve(i) for i in ids)

    def get_up_edges(self, node_id):
        return list(self._up.get(self._resolve(node_id), ()))

    def get_down_edges(self, node_id):
        return list(self._down.get(self._resolve(node_id), ()))

    def make_node(self, negative=False):
        # Real grgl::MutableGRG::makeNode always returns the plain,
        # non-negative slot id -- `negative` is not encoded in the return
        # value (see module docstring). Accepted here only for API parity.
        nid = self._next_node_id
        self._next_node_id += 1
        self._up[nid] = []
        self._down[nid] = []
        self._muts_by_node[nid] = []
        return nid

    def connect(self, parent, child):
        p, c = self._resolve(parent), self._resolve(child)
        self._up[c].append(p)
        self._down[p].append(c)

    def get_mutations_for_node(self, node_id):
        return list(self._muts_by_node.get(self._resolve(node_id), ()))

    def get_mutation_by_id(self, mut_id):
        return self._mutations[mut_id]

    def add_mutation(self, mutation, node_id):
        self._muts_by_node[self._resolve(node_id)].append(mutation.id)

    def remove_mutation(self, mut_id, node_id):
        self._muts_by_node[self._resolve(node_id)].remove(mut_id)

    def sort_mutations(self):
        for nid, ids in self._muts_by_node.items():
            ids.sort(key=lambda m: self._mutations[m].position)

    # ---- test-only helpers (not part of the pygrgl surface) ----

    def add_seed_mutation(self, node_id, position):
        mid = self._next_mut_id
        self._next_mut_id += 1
        self._mutations[mid] = _Mutation(mid, position)
        self._muts_by_node[node_id].append(mid)
        return mid

    def deepcopy(self):
        import copy
        return copy.deepcopy(self)

    def ancestor_mutation_positions(self, node_id):
        """Ground truth: positions of every mutation stored at node_id or
        any of its ancestors (BFS up via get_up_edges)."""
        seen = {node_id}
        queue = deque([node_id])
        positions = set()
        while queue:
            n = queue.popleft()
            for mid in self._muts_by_node.get(n, ()):
                positions.add(self._mutations[mid].position)
            for p in self._up.get(n, ()):
                if p not in seen:
                    seen.add(p)
                    queue.append(p)
        return positions


def _build_random_grg(rng, n_leaves=16, n_internal_layers=4, bp=200, mut_density=0.5):
    """A layered DAG (each layer's nodes connect up to ~2 nodes in the layer
    above, so nodes can have multiple parents / multiple children -- not a
    strict tree) with random mutations scattered at internal nodes and
    leaves alike."""
    g = FakeGRG((0, bp))
    layers = [[g.make_node() for _ in range(n_leaves)]]
    width = n_leaves
    for _ in range(n_internal_layers):
        width = max(1, width // 2)
        layer = [g.make_node() for _ in range(width)]
        prev = layers[-1]
        for i, child in enumerate(prev):
            targets = {layer[i % width], layer[(i + 1) % width]}
            for t in targets:
                g.connect(t, child)
        layers.append(layer)

    for layer in layers:
        for node in layer:
            if rng.random() < mut_density:
                for _ in range(rng.randint(1, 3)):
                    pos = rng.randint(0, bp - 1)
                    g.add_seed_mutation(node, pos)

    leaves = layers[0]
    g.set_samples(leaves)
    return g, leaves


def _random_segments(rng, leaves, bp, max_breakpoints=3):
    h1, h2 = rng.sample(leaves, 2)
    n_bp = rng.randint(0, max_breakpoints)
    cuts = sorted(rng.sample(range(1, bp), n_bp)) if n_bp else []
    parents = [h1, h2]
    segments = []
    start_parent = rng.randint(0, 1)
    prev = 0
    for i, c in enumerate(cuts):
        segments.append((parents[(start_parent + i) % 2], c))
        prev = c
    segments.append((parents[(start_parent + len(cuts)) % 2], bp))
    return segments


def _offspring_genotype(recomb, offspring_handle):
    """`recombine_multi` returns a small negative *handle* (an index into
    `recomb.NEGATIVE_NODE_IDS`, see `_register_offspring`), not the raw
    value `make_node(negative=True)` produced. `NEGATIVE_NODE_IDS` stores
    the real (non-negative) node id directly, so no further resolution
    is needed."""
    real_id = recomb.NEGATIVE_NODE_IDS[abs(offspring_handle) - 1]
    return recomb.grg.ancestor_mutation_positions(real_id)


def _expected_genotype(g_before, segments):
    expected = set()
    start = 0
    for parent_id, end in segments:
        if end > start:
            positions = g_before.ancestor_mutation_positions(parent_id)
            expected.update(p for p in positions if start <= p < end)
        start = end
    return expected


@pytest.mark.parametrize("seed", range(8))
def test_parallel_matches_sequential_and_ground_truth(seed):
    rng = random.Random(seed)
    g_template, leaves = _build_random_grg(rng, n_leaves=20, n_internal_layers=4, bp=200)

    all_segments = [_random_segments(rng, leaves, 200) for _ in range(15)]

    # Run sequential on its own copy.
    g_seq = g_template.deepcopy()
    recomb_seq = NonDuplicationRecombination(g_seq)
    seq_offspring = [recomb_seq.recombine_multi(segs) for segs in all_segments]
    recomb_seq.audit_check(raise_on_fail=True)

    # Run parallel on a fresh identical copy, single-worker (deterministic
    # thread interleaving isn't required for correctness, but keeping this
    # small avoids spurious flakiness from over-subscription in CI).
    g_par = g_template.deepcopy()
    recomb_par = ParallelNonDuplicationRecombination(g_par, max_workers=4)
    try:
        par_offspring = [recomb_par.recombine_multi(segs) for segs in all_segments]
        recomb_par.audit_check(raise_on_fail=True)
    finally:
        recomb_par.shutdown()

    for i, segs in enumerate(all_segments):
        expected = _expected_genotype(g_template, segs)
        seq_got = _offspring_genotype(recomb_seq, seq_offspring[i])
        par_got = _offspring_genotype(recomb_par, par_offspring[i])
        assert seq_got == expected, f"seed={seed} offspring={i}: sequential mismatch"
        assert par_got == expected, f"seed={seed} offspring={i}: parallel mismatch"


@pytest.mark.parametrize("seed", range(4))
def test_parallel_via_simulate_grg_recombination(seed):
    """End-to-end through the same driver used by the benchmarks."""
    rng = random.Random(1000 + seed)
    g_template, leaves = _build_random_grg(rng, n_leaves=24, n_internal_layers=4, bp=150)

    class Provider:
        num_offspring_per_couple = 2

        def get_breakpoints(self, bp_range, expected_crossovers=1.5):
            import numpy as np
            n = rng.randint(0, 3)
            if n == 0:
                return np.array([], dtype=int), 0
            low, high = bp_range
            bp = sorted(rng.sample(range(low + 1, high), n))
            return np.array(bp, dtype=int), n

    g = g_template.deepcopy()
    recomb = ParallelNonDuplicationRecombination(g, max_workers=4)
    try:
        offspring_ids, total_bp = simulate_grg_recombination(
            Provider(), recomb, g_template.bp_range, N=g_template.bp_range[1]
        )
        recomb.audit_check(raise_on_fail=True)
    finally:
        recomb.shutdown()

    assert len(offspring_ids) > 0
    assert set(g.get_sample_nodes()) == set(g._resolve(i) for i in offspring_ids)
