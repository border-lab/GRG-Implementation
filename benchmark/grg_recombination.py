"""
grg_recombination.py

Optimized core module containing the Non-duplication GRG recombination algorithm.
"""

import bisect
import time
from collections import deque

import numpy as np


def _summary(xs):
    """min/mean/p50/p95/p99/max for a list of numbers (or {'count': 0} if empty)."""
    if not xs:
        return {"count": 0}
    xs_sorted = sorted(xs)
    m = len(xs_sorted)
    return {
        "count": m,
        "min": xs_sorted[0],
        "max": xs_sorted[-1],
        "mean": sum(xs_sorted) / m,
        "p50": xs_sorted[m // 2],
        "p95": xs_sorted[min(m - 1, int(m * 0.95))],
        "p99": xs_sorted[min(m - 1, int(m * 0.99))],
    }


def compute_grg_structural_stats(grg):
    """One-shot graph topology metrics; safe to run before recombination.
    O(V + E) plus sorts for distribution summaries."""
    n = grg.num_nodes
    samples = [i for i in range(n) if grg.is_sample(i)]
    roots = list(grg.get_root_nodes())

    # Topological depth via BFS down from roots.
    depth = {r: 0 for r in roots}
    queue = deque(roots)
    while queue:
        node = queue.popleft()
        d = depth[node]
        for child in grg.get_down_edges(node):
            cd = depth.get(child, -1)
            if cd < d + 1:
                depth[child] = d + 1
                queue.append(child)
    sample_depths = [depth.get(s, 0) for s in samples]

    mut_counts = [len(grg.get_mutations_for_node(i)) for i in range(n)]
    up_fanouts = [len(grg.get_up_edges(i)) for i in range(n)]
    down_fanouts = [len(grg.get_down_edges(i)) for i in range(n)]

    return {
        "num_nodes": n,
        "num_edges": grg.num_edges,
        "num_samples": len(samples),
        "num_mutations": grg.num_mutations,
        "num_roots": len(roots),
        "internal_per_sample": (n - len(samples)) / max(1, len(samples)),
        "genome_length": grg.bp_range[1],
        "mutations_per_node": _summary(mut_counts),
        "up_fanout": _summary(up_fanouts),
        "down_fanout": _summary(down_fanouts),
        "sample_depth": _summary(sample_depths),
    }

def recombination_intervals(h1, h2, bp, N):
    """
    Returns list of segments: [(source_parent_id, end_coord), ...]
    """
    # bp = self.get_breakpoints(N)
    start = np.random.binomial(1, 0.5, 1)[0]
    parents = [h1, h2]

    segments = []
    for i, K in enumerate(bp):
        segments.append((parents[(start + i) % 2], K))
    segments.append((parents[(start + len(bp)) % 2], N))
    return segments


class NonDuplicationRecombination:
    """
    Non-duplication GRG recombination algorithm.

    Optimizations vs. the recursive version:
    - Iterative DFS for `_recurse_attach` and `_get_node_and_ancestor_span`
      (avoids Python's recursion limit on deep GRGs).
    - O(1) offspring reverse lookup via a dict alongside NEGATIVE_NODE_IDS.
    - Per-recombiner cache for `get_up_edges()`.
    - Generation-counter arrays replace per-call Python sets for
      `visited` (per DFS) and `connected` (per offspring).
    - Optional deferred sample updates via `defer_sample_updates`.
    - Hot-loop micro-optimizations: cache-hit fast paths for mutation
      range and ancestral coverage are inlined into `_recurse_attach`,
      with class attributes hoisted into local variables (CPython
      optimizes LOAD_FAST faster than LOAD_ATTR).
    - `_extract_bubble` no longer invalidates the caches of node_id's
      parents; bubble extraction does not affect their span or ancestral
      coverage, so invalidation just thrashes high-traffic cache entries.
    - Per-parent interval batching: `recombine_multi` groups all segments
      contributed by the same parent into a single traversal whose query
      is the union of those intervals. This avoids re-bubbling the same
      ancestor when one parent contributes mutations across multiple
      disjoint segments of the offspring's genome.
    """

    debug_mode = False

    def __init__(self, grg, instrument=False):
        self.grg = grg
        self.genome_length = grg.bp_range[1]
        self.original_bp_range = grg.bp_range
        self.NEGATIVE_NODE_IDS = []
        self._negative_node_index = {}
        self._modified_nodes = set()
        self._pending_bubbles = []
        self._pending_sample_removals = set()

        self.defer_sample_updates = False

        self._mutation_cache = {}
        self._pos_cache = {}
        self._up_edges_cache = {}

        self.span_cache = [False] * self.grg.num_nodes
        self.anc_cov_cache = [False] * self.grg.num_nodes

        self._visited_gen = [0] * self.grg.num_nodes
        self._connected_gen = [0] * self.grg.num_nodes
        self._gen_visited = 0
        self._gen_connected = 0

        # ----- Instrumentation -----
        # When True, every public recombine call accumulates phase / C++-call
        # timings and per-offspring + per-bubble distributions into self.stats.
        # Overhead per moved mutation is ~150 ns (3 perf_counter calls) -- well
        # below 0.1% on the workloads we care about.
        self.instrument = instrument
        self.stats = self._fresh_stats()

    @staticmethod
    def _fresh_stats():
        return {
            # Phase-level wallclock totals (seconds)
            "recurse_attach_time": 0.0,
            "apply_bubbles_time": 0.0,
            "sync_to_grg_time": 0.0,
            "clear_caches_time": 0.0,
            "flush_samples_time": 0.0,
            # C++-call wallclock totals (seconds)
            "get_mutation_by_id_time": 0.0,
            "add_mutation_time": 0.0,
            "remove_mutation_time": 0.0,
            "make_node_time": 0.0,
            "connect_time": 0.0,
            # C++-call counts (per_call_cost = time / count)
            "get_mutation_by_id_calls": 0,
            "add_mutation_calls": 0,
            "remove_mutation_calls": 0,
            "make_node_calls": 0,
            "connect_calls": 0,
            # Algorithmic counters
            "offspring_count": 0,
            "recurse_attach_calls": 0,
            "segments_processed": 0,
            "visits_total": 0,
            "bubbles_created": 0,
            "mutations_moved": 0,
        }

    def reset_stats(self):
        """Clear instrumentation accumulators (e.g. between generations)."""
        self.stats = self._fresh_stats()

    # ------------------------------------------------------------------
    # Per-node array growth
    # ------------------------------------------------------------------

    def _grow_node_arrays(self, node_id):
        target = node_id + 1
        n = len(self.span_cache)
        if n < target:
            pad = target - n
            self.span_cache.extend([False] * pad)
            self.anc_cov_cache.extend([False] * pad)
        n = len(self._visited_gen)
        if n < target:
            pad = target - n
            self._visited_gen.extend([0] * pad)
            self._connected_gen.extend([0] * pad)

    def _sync_to_grg(self):
        n = self.grg.num_nodes
        if n > 0:
            self._grow_node_arrays(n - 1)

    # ------------------------------------------------------------------
    # Cached graph access
    # ------------------------------------------------------------------

    def _get_up_edges_cached(self, node_id):
        cached = self._up_edges_cache.get(node_id)
        if cached is None:
            cached = list(self.grg.get_up_edges(node_id))
            self._up_edges_cache[node_id] = cached
        return cached

    # ------------------------------------------------------------------
    # Mutation caches
    # ------------------------------------------------------------------

    def _get_node_mutations(self, node_id):
        if node_id not in self._mutation_cache:
            mut_ids = self.grg.get_mutations_for_node(node_id)
            mutations = []
            for mut_id in mut_ids:
                mut = self.grg.get_mutation_by_id(mut_id)
                mutations.append((mut_id, mut.position))
            combined = sorted(mutations, key=lambda x: x[1])
            self._mutation_cache[node_id] = combined
            self._pos_cache[node_id] = [m[1] for m in combined]
        return self._mutation_cache[node_id]

    def _get_mutation_range(self, node_id, L, R):
        self._get_node_mutations(node_id)
        positions = self._pos_cache[node_id]
        if not positions:
            return 0, 0
        return bisect.bisect_left(positions, L), bisect.bisect_left(positions, R)

    # ------------------------------------------------------------------
    # Interval-union helpers
    # ------------------------------------------------------------------
    # `intervals` is always a sorted tuple of disjoint, coalesced,
    # half-open (L, R) pairs representing the union of genomic regions
    # the offspring still needs to inherit on the current traversal.

    @staticmethod
    def _intervals_disjoint_from(intervals, lo, hi):
        for L, R in intervals:
            if L >= hi:
                return True
            if R > lo:
                return False
        return True

    @staticmethod
    def _intervals_contain(intervals, lo, hi):
        # True iff some single piece (L_i, R_i) of the union has
        # L_i <= lo and hi <= R_i. Because the pieces are disjoint and
        # coalesced, no half-open [lo, hi) can be covered by spanning
        # multiple pieces -- coalescing collapses adjacency.
        for L, R in intervals:
            if L > lo:
                return False
            if R >= hi:
                return True
        return False

    @staticmethod
    def _clip_intervals(intervals, lo, hi):
        result = []
        for L, R in intervals:
            if R <= lo:
                continue
            if L >= hi:
                break
            nL = L if L > lo else lo
            nR = R if R < hi else hi
            result.append((nL, nR))
        return tuple(result)

    # ------------------------------------------------------------------
    # Iterative span / ancestral coverage
    # ------------------------------------------------------------------

    def _get_node_and_ancestor_span(self, node_id):
        if self.span_cache[node_id] is not False:
            return self.span_cache[node_id]

        stack = [(node_id, False)]
        scheduled = {node_id}

        while stack:
            nid, processed = stack.pop()

            if processed:
                min_p, max_p = float('inf'), float('-inf')
                node_muts = self._get_node_mutations(nid)
                if node_muts:
                    min_p = node_muts[0][1]
                    max_p = node_muts[-1][1]
                for parent in self._get_up_edges_cached(nid):
                    anc = self.span_cache[parent]
                    if anc:
                        if anc[0] < min_p: min_p = anc[0]
                        if anc[1] > max_p: max_p = anc[1]
                self.span_cache[nid] = None if min_p == float('inf') else (min_p, max_p)
                scheduled.discard(nid)
            else:
                stack.append((nid, True))
                for parent in self._get_up_edges_cached(nid):
                    if self.span_cache[parent] is False and parent not in scheduled:
                        stack.append((parent, False))
                        scheduled.add(parent)

        return self.span_cache[node_id]

    def _get_ancestral_coverage(self, node_id):
        if self.anc_cov_cache[node_id] is not False:
            return self.anc_cov_cache[node_id]

        parents = self._get_up_edges_cached(node_id)
        if not parents:
            self.anc_cov_cache[node_id] = None
            return None

        min_pos, max_pos = float('inf'), float('-inf')
        for parent in parents:
            p_span = self._get_node_and_ancestor_span(parent)
            if p_span:
                if p_span[0] < min_pos: min_pos = p_span[0]
                if p_span[1] > max_pos: max_pos = p_span[1]

        result = None if min_pos == float('inf') else (min_pos, max_pos + 1)
        self.anc_cov_cache[node_id] = result
        return result

    # ------------------------------------------------------------------
    # Bubble extraction
    # ------------------------------------------------------------------

    def _extract_bubble(self, node_id, relevant_mut_ids, offspring_id, intervals):
        instrument = self.instrument
        if instrument:
            t = time.perf_counter()
        bubble_id = self.grg.make_node()
        if instrument:
            self.stats["make_node_time"] += time.perf_counter() - t
            self.stats["make_node_calls"] += 1

        self._grow_node_arrays(bubble_id)

        if instrument:
            t = time.perf_counter()
        self.grg.connect(bubble_id, node_id)
        self.grg.connect(bubble_id, -offspring_id)
        if instrument:
            self.stats["connect_time"] += time.perf_counter() - t
            self.stats["connect_calls"] += 2

        # node_id just gained a new up-edge; drop its cached up-edges list.
        self._up_edges_cache.pop(node_id, None)

        self._pending_bubbles.append({
            'node_id': node_id,
            'bubble_id': bubble_id,
            'relevant_mut_ids': relevant_mut_ids,
        })

        # Track only nodes whose own caches actually become stale:
        # node_id loses mutations (mutation/pos/anc_cov are stale) and the
        # new bubble_id needs a clean slot. node_id's parents are NOT
        # invalidated -- their span and anc_cov depend only on themselves
        # and their own ancestors, neither of which changes here.
        self._modified_nodes.add(node_id)
        self._modified_nodes.add(bubble_id)

        return bubble_id

    # ------------------------------------------------------------------
    # Iterative recurse-attach (hot path)
    # ------------------------------------------------------------------

    def _recurse_attach(self, root_id, offspring_id, intervals):
        """`intervals`: sorted, coalesced tuple of disjoint half-open (L, R) pairs."""
        if not intervals:
            return

        self._gen_visited += 1
        gen_v = self._gen_visited
        gen_c = self._gen_connected

        # Hoist attributes / globals into locals for the inner loop.
        # CPython's LOAD_FAST is faster than LOAD_ATTR; over hundreds of
        # thousands of iterations these adds up.
        visited_gen = self._visited_gen
        connected_gen = self._connected_gen
        pos_cache = self._pos_cache
        mutation_cache = self._mutation_cache
        anc_cov_cache = self.anc_cov_cache
        bisect_left = bisect.bisect_left
        get_node_mutations = self._get_node_mutations
        get_ancestral_coverage = self._get_ancestral_coverage
        get_up_edges_cached = self._get_up_edges_cached
        extract_bubble = self._extract_bubble
        clip_intervals = NonDuplicationRecombination._clip_intervals
        grg_connect = self.grg.connect
        pending_sample_removals_add = self._pending_sample_removals.add

        neg_offspring = -offspring_id
        stack = [(root_id, intervals)]
        stack_append = stack.append
        stack_pop = stack.pop

        while stack:
            node_id, ivs = stack_pop()

            if not ivs:
                continue
            if visited_gen[node_id] == gen_v:
                continue
            visited_gen[node_id] = gen_v

            # Multi-bisect over each piece of the union; collect both the
            # count of relevant mutations and the (left, right) slice
            # bounds (only used if we end up extracting a bubble).
            positions = pos_cache.get(node_id)
            if positions is None:
                get_node_mutations(node_id)              # populates both caches
                positions = pos_cache[node_id]

            num_all = len(positions)
            num_rel = 0
            slice_bounds = None
            if num_all:
                slice_bounds = []
                for L, R in ivs:
                    left = bisect_left(positions, L)
                    right = bisect_left(positions, R)
                    if right > left:
                        slice_bounds.append((left, right))
                        num_rel += right - left

            has_all_relevant = num_all > 0 and num_rel == num_all
            has_no_relevant = num_rel == 0
            has_partial_relevant = num_rel > 0 and num_rel < num_all

            # Inlined _get_ancestral_coverage cache-hit path.
            Iu = anc_cov_cache[node_id]
            if Iu is False:
                Iu = get_ancestral_coverage(node_id)

            if Iu is None:
                if has_all_relevant and connected_gen[node_id] != gen_c:
                    grg_connect(node_id, neg_offspring)
                    connected_gen[node_id] = gen_c
                    pending_sample_removals_add(node_id)
                if has_partial_relevant:
                    muts = mutation_cache[node_id]
                    rel_mut_ids = []
                    for l, r in slice_bounds:
                        rel_mut_ids.extend(m[0] for m in muts[l:r])
                    bubble_id = extract_bubble(node_id, rel_mut_ids, offspring_id, ivs)
                    connected_gen[bubble_id] = gen_c
                continue

            Iu0 = Iu[0]
            Iu1 = Iu[1]

            # ancestral_disjoint: no piece of the union overlaps [Iu0, Iu1).
            ancestral_disjoint = True
            for L, R in ivs:
                if L >= Iu1:
                    break
                if R > Iu0:
                    ancestral_disjoint = False
                    break

            # full_coverage: some single piece of the union fully contains
            # [Iu0, Iu1). Coalesced disjoint pieces can't cover [Iu0, Iu1)
            # by spanning, so a single-piece check is sufficient.
            full_coverage = False
            for L, R in ivs:
                if L > Iu0:
                    break
                if R >= Iu1:
                    full_coverage = True
                    break

            if has_all_relevant and full_coverage:
                if connected_gen[node_id] != gen_c:
                    grg_connect(node_id, neg_offspring)
                    connected_gen[node_id] = gen_c
                    pending_sample_removals_add(node_id)
                continue

            if has_no_relevant and ancestral_disjoint:
                continue

            if has_no_relevant:
                clipped = clip_intervals(ivs, Iu0, Iu1)
                if not clipped:
                    continue
                for parent in reversed(get_up_edges_cached(node_id)):
                    stack_append((parent, clipped))
                continue

            # Partial / all relevant, not full coverage -> bubble + maybe recurse.
            muts = mutation_cache[node_id]
            rel_mut_ids = []
            for l, r in slice_bounds:
                rel_mut_ids.extend(m[0] for m in muts[l:r])
            bubble_id = extract_bubble(node_id, rel_mut_ids, offspring_id, ivs)
            connected_gen[bubble_id] = gen_c

            if not ancestral_disjoint:
                clipped = clip_intervals(ivs, Iu0, Iu1)
                if not clipped:
                    continue
                for parent in reversed(get_up_edges_cached(node_id)):
                    stack_append((parent, clipped))

    # ------------------------------------------------------------------
    # Apply deferred work, evict caches
    # ------------------------------------------------------------------

    def _apply_pending_bubbles(self):
        instrument = self.instrument
        if instrument:
            t_apply_start = time.perf_counter()
            t_get = 0.0
            t_add = 0.0
            t_rem = 0.0
            n_bubbles = len(self._pending_bubbles)
            n_muts_total = 0

        for bubble_op in self._pending_bubbles:
            node_id = bubble_op['node_id']
            bubble_id = bubble_op['bubble_id']
            muts = bubble_op['relevant_mut_ids']
            if instrument:
                n_muts_total += len(muts)
            for mut_id in muts:
                if instrument:
                    t = time.perf_counter()
                mut = self.grg.get_mutation_by_id(mut_id)
                if instrument:
                    t_get += time.perf_counter() - t
                    t = time.perf_counter()
                self.grg.add_mutation(mut, bubble_id)
                if instrument:
                    t_add += time.perf_counter() - t
                    t = time.perf_counter()
                self.grg.remove_mutation(mut_id, node_id)
                if instrument:
                    t_rem += time.perf_counter() - t
        self._pending_bubbles.clear()

        if not self.defer_sample_updates:
            self.flush_sample_updates()

        if instrument:
            elapsed = time.perf_counter() - t_apply_start
            self.stats["apply_bubbles_time"] += elapsed
            self.stats["get_mutation_by_id_time"] += t_get
            self.stats["add_mutation_time"] += t_add
            self.stats["remove_mutation_time"] += t_rem
            self.stats["get_mutation_by_id_calls"] += n_muts_total
            self.stats["add_mutation_calls"] += n_muts_total
            self.stats["remove_mutation_calls"] += n_muts_total
            self.stats["bubbles_created"] += n_bubbles
            self.stats["mutations_moved"] += n_muts_total

    def flush_sample_updates(self):
        if self._pending_sample_removals:
            instrument = self.instrument
            if instrument:
                t = time.perf_counter()
            current = set(self.grg.get_sample_nodes())
            current.difference_update(self._pending_sample_removals)
            self.grg.set_samples(list(current))
            self._pending_sample_removals.clear()
            if instrument:
                self.stats["flush_samples_time"] += time.perf_counter() - t

    def _clear_modified_caches(self):
        if self.instrument:
            t = time.perf_counter()
        for node_id in self._modified_nodes:
            self._mutation_cache.pop(node_id, None)
            self._pos_cache.pop(node_id, None)
            if node_id < len(self.span_cache):
                self.span_cache[node_id] = False
                self.anc_cov_cache[node_id] = False
        self._modified_nodes.clear()
        if self.instrument:
            self.stats["clear_caches_time"] += time.perf_counter() - t

    # ------------------------------------------------------------------
    # Public entry points
    # ------------------------------------------------------------------

    def _register_offspring(self, offspring_id):
        idx = self._negative_node_index.get(offspring_id)
        if idx is None:
            idx = len(self.NEGATIVE_NODE_IDS)
            self._negative_node_index[offspring_id] = idx
            self.NEGATIVE_NODE_IDS.append(offspring_id)
        return -(idx + 1)

    def recombine(self, haplotype_A, haplotype_B, breakpoint):
        instrument = self.instrument
        self._pending_bubbles.clear()

        if instrument:
            t = time.perf_counter()
            gen_v_before = self._gen_visited
        self._sync_to_grg()
        if instrument:
            self.stats["sync_to_grg_time"] += time.perf_counter() - t
            t = time.perf_counter()
        offspring_id = self.grg.make_node(negative=True)
        if instrument:
            self.stats["make_node_time"] += time.perf_counter() - t
            self.stats["make_node_calls"] += 1
        self._grow_node_arrays(self.grg.num_nodes - 1)
        self._gen_connected += 1

        if instrument:
            recurse_t = 0.0
            ts = time.perf_counter()
        self._recurse_attach(haplotype_A, offspring_id, ((0, breakpoint),))
        if instrument:
            recurse_t += time.perf_counter() - ts
            self.stats["recurse_attach_calls"] += 1
            self.stats["segments_processed"] += 1
            ts = time.perf_counter()
        self._recurse_attach(haplotype_B, offspring_id, ((breakpoint, self.genome_length),))
        if instrument:
            recurse_t += time.perf_counter() - ts
            self.stats["recurse_attach_calls"] += 1
            self.stats["segments_processed"] += 1
            self.stats["recurse_attach_time"] += recurse_t

        self._apply_pending_bubbles()
        self._clear_modified_caches()

        if instrument:
            gen_v_after = self._gen_visited
            gen_set = set(range(gen_v_before + 1, gen_v_after + 1))
            v = sum(1 for x in self._visited_gen if x in gen_set)
            self.stats["visits_total"] += v
            self.stats["offspring_count"] += 1

        return self._register_offspring(offspring_id)

    def recombine_multi(self, segments):
        instrument = self.instrument
        self._pending_bubbles.clear()

        if instrument:
            t = time.perf_counter()
            gen_v_before = self._gen_visited
        self._sync_to_grg()
        if instrument:
            self.stats["sync_to_grg_time"] += time.perf_counter() - t
            t = time.perf_counter()
        offspring_id = self.grg.make_node(negative=True)
        if instrument:
            self.stats["make_node_time"] += time.perf_counter() - t
            self.stats["make_node_calls"] += 1
        self._grow_node_arrays(self.grg.num_nodes - 1)
        self._gen_connected += 1

        # Bucket segments by parent so each parent gets one traversal
        # over the union of its contributed intervals. Avoids re-bubbling
        # the same ancestor when one parent contributes mutations to two
        # separated regions of the offspring's genome.
        by_parent = {}
        start = 0
        seg_count = 0
        for parent_id, end in segments:
            if end > start:
                by_parent.setdefault(parent_id, []).append((start, end))
                seg_count += 1
            start = end

        if instrument:
            recurse_t = 0.0
        for parent_id, ivs in by_parent.items():
            ivs.sort()
            # Coalesce adjacent / overlapping intervals so the coverage
            # checks in `_recurse_attach` (which assume coalesced input)
            # remain correct. `recombination_intervals` produces strictly
            # alternating segments so adjacency is rare, but selfing or
            # odd breakpoint patterns can produce it.
            merged = []
            for L, R in ivs:
                if merged and L <= merged[-1][1]:
                    prev_L, prev_R = merged[-1]
                    merged[-1] = (prev_L, R if R > prev_R else prev_R)
                else:
                    merged.append((L, R))
            if self.debug_mode:
                print(f"BREAK parent={parent_id} intervals={merged}")
            if instrument:
                ts = time.perf_counter()
            self._recurse_attach(parent_id, offspring_id, tuple(merged))
            if instrument:
                recurse_t += time.perf_counter() - ts
                self.stats["recurse_attach_calls"] += 1
        if instrument:
            self.stats["recurse_attach_time"] += recurse_t
            self.stats["segments_processed"] += seg_count

        self._apply_pending_bubbles()
        self._clear_modified_caches()

        if instrument:
            gen_v_after = self._gen_visited
            gen_set = set(range(gen_v_before + 1, gen_v_after + 1))
            v = sum(1 for x in self._visited_gen if x in gen_set)
            self.stats["visits_total"] += v
            self.stats["offspring_count"] += 1

        return self._register_offspring(offspring_id)


class NonDuplicationRecombinationLegacy:
    """
    Non-duplication GRG recombination algorithm.

    Optimizations vs. the recursive version:
    - Iterative DFS for `_recurse_attach` and `_get_node_and_ancestor_span`
      (avoids Python's recursion limit on deep GRGs).
    - O(1) offspring reverse lookup via a dict alongside NEGATIVE_NODE_IDS.
    - Per-recombiner cache for `get_up_edges()`.
    - Generation-counter arrays replace per-call Python sets for
      `visited` (per DFS) and `connected` (per offspring).
    - Optional deferred sample updates via `defer_sample_updates`.
    - Hot-loop micro-optimizations: cache-hit fast paths for mutation
      range and ancestral coverage are inlined into `_recurse_attach`,
      with class attributes hoisted into local variables (CPython
      optimizes LOAD_FAST faster than LOAD_ATTR).
    - `_extract_bubble` no longer invalidates the caches of node_id's
      parents; bubble extraction does not affect their span or ancestral
      coverage, so invalidation just thrashes high-traffic cache entries.
    """

    debug_mode = False

    def __init__(self, grg, instrument=False):
        self.grg = grg
        self.genome_length = grg.bp_range[1]
        self.original_bp_range = grg.bp_range
        self.NEGATIVE_NODE_IDS = []
        self._negative_node_index = {}
        self._modified_nodes = set()
        self._pending_bubbles = []
        self._pending_sample_removals = set()

        self.defer_sample_updates = False

        self._mutation_cache = {}
        self._pos_cache = {}
        self._up_edges_cache = {}

        self.span_cache = [False] * self.grg.num_nodes
        self.anc_cov_cache = [False] * self.grg.num_nodes

        self._visited_gen = [0] * self.grg.num_nodes
        self._connected_gen = [0] * self.grg.num_nodes
        self._gen_visited = 0
        self._gen_connected = 0

        # ----- Instrumentation -----
        # When True, every public recombine call accumulates phase / C++-call
        # timings and per-offspring + per-bubble distributions into self.stats.
        # Overhead per moved mutation is ~150 ns (3 perf_counter calls) -- well
        # below 0.1% on the workloads we care about.
        self.instrument = instrument
        self.stats = self._fresh_stats()

    @staticmethod
    def _fresh_stats():
        return {
            # Phase-level wallclock totals (seconds)
            "recurse_attach_time": 0.0,
            "apply_bubbles_time": 0.0,
            "sync_to_grg_time": 0.0,
            "clear_caches_time": 0.0,
            "flush_samples_time": 0.0,
            # C++-call wallclock totals (seconds)
            "get_mutation_by_id_time": 0.0,
            "add_mutation_time": 0.0,
            "remove_mutation_time": 0.0,
            "make_node_time": 0.0,
            "connect_time": 0.0,
            # C++-call counts (per_call_cost = time / count)
            "get_mutation_by_id_calls": 0,
            "add_mutation_calls": 0,
            "remove_mutation_calls": 0,
            "make_node_calls": 0,
            "connect_calls": 0,
            # Algorithmic counters
            "offspring_count": 0,
            "recurse_attach_calls": 0,
            "segments_processed": 0,
            "visits_total": 0,
            "bubbles_created": 0,
            "mutations_moved": 0,
        }

    def reset_stats(self):
        """Clear instrumentation accumulators (e.g. between generations)."""
        self.stats = self._fresh_stats()

    # ------------------------------------------------------------------
    # Per-node array growth
    # ------------------------------------------------------------------

    def _grow_node_arrays(self, node_id):
        target = node_id + 1
        n = len(self.span_cache)
        if n < target:
            pad = target - n
            self.span_cache.extend([False] * pad)
            self.anc_cov_cache.extend([False] * pad)
        n = len(self._visited_gen)
        if n < target:
            pad = target - n
            self._visited_gen.extend([0] * pad)
            self._connected_gen.extend([0] * pad)

    def _sync_to_grg(self):
        n = self.grg.num_nodes
        if n > 0:
            self._grow_node_arrays(n - 1)

    # ------------------------------------------------------------------
    # Cached graph access
    # ------------------------------------------------------------------

    def _get_up_edges_cached(self, node_id):
        cached = self._up_edges_cache.get(node_id)
        if cached is None:
            cached = list(self.grg.get_up_edges(node_id))
            self._up_edges_cache[node_id] = cached
        return cached

    # ------------------------------------------------------------------
    # Mutation caches
    # ------------------------------------------------------------------

    def _get_node_mutations(self, node_id):
        if node_id not in self._mutation_cache:
            mut_ids = self.grg.get_mutations_for_node(node_id)
            mutations = []
            for mut_id in mut_ids:
                mut = self.grg.get_mutation_by_id(mut_id)
                mutations.append((mut_id, mut.position))
            combined = sorted(mutations, key=lambda x: x[1])
            self._mutation_cache[node_id] = combined
            self._pos_cache[node_id] = [m[1] for m in combined]
        return self._mutation_cache[node_id]

    def _get_mutation_range(self, node_id, L, R):
        self._get_node_mutations(node_id)
        positions = self._pos_cache[node_id]
        if not positions:
            return 0, 0
        return bisect.bisect_left(positions, L), bisect.bisect_left(positions, R)

    # ------------------------------------------------------------------
    # Iterative span / ancestral coverage
    # ------------------------------------------------------------------

    def _get_node_and_ancestor_span(self, node_id):
        if self.span_cache[node_id] is not False:
            return self.span_cache[node_id]

        stack = [(node_id, False)]
        scheduled = {node_id}

        while stack:
            nid, processed = stack.pop()

            if processed:
                min_p, max_p = float('inf'), float('-inf')
                node_muts = self._get_node_mutations(nid)
                if node_muts:
                    min_p = node_muts[0][1]
                    max_p = node_muts[-1][1]
                for parent in self._get_up_edges_cached(nid):
                    anc = self.span_cache[parent]
                    if anc:
                        if anc[0] < min_p: min_p = anc[0]
                        if anc[1] > max_p: max_p = anc[1]
                self.span_cache[nid] = None if min_p == float('inf') else (min_p, max_p)
                scheduled.discard(nid)
            else:
                stack.append((nid, True))
                for parent in self._get_up_edges_cached(nid):
                    if self.span_cache[parent] is False and parent not in scheduled:
                        stack.append((parent, False))
                        scheduled.add(parent)

        return self.span_cache[node_id]

    def _get_ancestral_coverage(self, node_id):
        if self.anc_cov_cache[node_id] is not False:
            return self.anc_cov_cache[node_id]

        parents = self._get_up_edges_cached(node_id)
        if not parents:
            self.anc_cov_cache[node_id] = None
            return None

        min_pos, max_pos = float('inf'), float('-inf')
        for parent in parents:
            p_span = self._get_node_and_ancestor_span(parent)
            if p_span:
                if p_span[0] < min_pos: min_pos = p_span[0]
                if p_span[1] > max_pos: max_pos = p_span[1]

        result = None if min_pos == float('inf') else (min_pos, max_pos + 1)
        self.anc_cov_cache[node_id] = result
        return result

    # ------------------------------------------------------------------
    # Bubble extraction
    # ------------------------------------------------------------------

    def _extract_bubble(self, node_id, relevant_mut_ids, offspring_id, interval):
        instrument = self.instrument
        if instrument:
            t = time.perf_counter()
        bubble_id = self.grg.make_node()
        if instrument:
            self.stats["make_node_time"] += time.perf_counter() - t
            self.stats["make_node_calls"] += 1

        self._grow_node_arrays(bubble_id)

        if instrument:
            t = time.perf_counter()
        self.grg.connect(bubble_id, node_id)
        self.grg.connect(bubble_id, -offspring_id)
        if instrument:
            self.stats["connect_time"] += time.perf_counter() - t
            self.stats["connect_calls"] += 2

        # node_id just gained a new up-edge; drop its cached up-edges list.
        self._up_edges_cache.pop(node_id, None)

        self._pending_bubbles.append({
            'node_id': node_id,
            'bubble_id': bubble_id,
            'relevant_mut_ids': relevant_mut_ids,
        })

        # Track only nodes whose own caches actually become stale:
        # node_id loses mutations (mutation/pos/anc_cov are stale) and the
        # new bubble_id needs a clean slot. node_id's parents are NOT
        # invalidated -- their span and anc_cov depend only on themselves
        # and their own ancestors, neither of which changes here.
        self._modified_nodes.add(node_id)
        self._modified_nodes.add(bubble_id)

        return bubble_id

    # ------------------------------------------------------------------
    # Iterative recurse-attach (hot path)
    # ------------------------------------------------------------------

    def _recurse_attach(self, root_id, offspring_id, L0, R0):
        self._gen_visited += 1
        gen_v = self._gen_visited
        gen_c = self._gen_connected

        # Hoist attributes / globals into locals for the inner loop.
        # CPython's LOAD_FAST is faster than LOAD_ATTR; over hundreds of
        # thousands of iterations these adds up.
        visited_gen = self._visited_gen
        connected_gen = self._connected_gen
        pos_cache = self._pos_cache
        mutation_cache = self._mutation_cache
        anc_cov_cache = self.anc_cov_cache
        bisect_left = bisect.bisect_left
        get_node_mutations = self._get_node_mutations
        get_ancestral_coverage = self._get_ancestral_coverage
        get_up_edges_cached = self._get_up_edges_cached
        extract_bubble = self._extract_bubble
        grg_connect = self.grg.connect
        pending_sample_removals_add = self._pending_sample_removals.add

        neg_offspring = -offspring_id
        stack = [(root_id, L0, R0)]
        stack_append = stack.append
        stack_pop = stack.pop

        while stack:
            node_id, L, R = stack_pop()

            if L >= R:
                continue
            if visited_gen[node_id] == gen_v:
                continue
            visited_gen[node_id] = gen_v

            # Inlined _get_mutation_range + _get_node_mutations cache-hit path.
            positions = pos_cache.get(node_id)
            if positions is None:
                get_node_mutations(node_id)              # populates both caches
                positions = pos_cache[node_id]

            if positions:
                left = bisect_left(positions, L)
                right = bisect_left(positions, R)
            else:
                left = right = 0

            num_rel = right - left
            num_all = len(positions)

            has_all_relevant = num_all > 0 and num_rel == num_all
            has_no_relevant = num_rel == 0
            has_partial_relevant = num_rel > 0 and num_rel < num_all

            # Inlined _get_ancestral_coverage cache-hit path.
            Iu = anc_cov_cache[node_id]
            if Iu is False:
                Iu = get_ancestral_coverage(node_id)

            if Iu is None:
                if has_all_relevant and connected_gen[node_id] != gen_c:
                    grg_connect(node_id, neg_offspring)
                    connected_gen[node_id] = gen_c
                    pending_sample_removals_add(node_id)
                if has_partial_relevant:
                    rel_mut_ids = [m[0] for m in mutation_cache[node_id][left:right]]
                    bubble_id = extract_bubble(node_id, rel_mut_ids, offspring_id, (L, R))
                    connected_gen[bubble_id] = gen_c
                continue

            Iu0 = Iu[0]
            Iu1 = Iu[1]
            ancestral_disjoint = R <= Iu0 or L >= Iu1
            full_coverage = Iu0 >= L and Iu1 <= R

            if has_all_relevant and full_coverage:
                if connected_gen[node_id] != gen_c:
                    grg_connect(node_id, neg_offspring)
                    connected_gen[node_id] = gen_c
                    pending_sample_removals_add(node_id)
                continue

            if has_no_relevant and ancestral_disjoint:
                continue

            if has_no_relevant:
                # Inline max/min — Python's builtins are surprisingly slow
                # at this call rate (showed up at ~50ms in the profile).
                newL = L if L > Iu0 else Iu0
                newR = R if R < Iu1 else Iu1
                if newL >= newR:
                    continue
                for parent in reversed(get_up_edges_cached(node_id)):
                    stack_append((parent, newL, newR))
                continue

            # Partial / all relevant, not full coverage -> bubble + maybe recurse.
            rel_mut_ids = [m[0] for m in mutation_cache[node_id][left:right]]
            bubble_id = extract_bubble(node_id, rel_mut_ids, offspring_id, (L, R))
            connected_gen[bubble_id] = gen_c

            if not ancestral_disjoint:
                newL = L if L > Iu0 else Iu0
                newR = R if R < Iu1 else Iu1
                if newL >= newR:
                    continue
                for parent in reversed(get_up_edges_cached(node_id)):
                    stack_append((parent, newL, newR))

    # ------------------------------------------------------------------
    # Apply deferred work, evict caches
    # ------------------------------------------------------------------

    def _apply_pending_bubbles(self):
        instrument = self.instrument
        if instrument:
            t_apply_start = time.perf_counter()
            t_get = 0.0
            t_add = 0.0
            t_rem = 0.0
            n_bubbles = len(self._pending_bubbles)
            n_muts_total = 0

        for bubble_op in self._pending_bubbles:
            node_id = bubble_op['node_id']
            bubble_id = bubble_op['bubble_id']
            muts = bubble_op['relevant_mut_ids']
            if instrument:
                n_muts_total += len(muts)
            for mut_id in muts:
                if instrument:
                    t = time.perf_counter()
                mut = self.grg.get_mutation_by_id(mut_id)
                if instrument:
                    t_get += time.perf_counter() - t
                    t = time.perf_counter()
                self.grg.add_mutation(mut, bubble_id)
                if instrument:
                    t_add += time.perf_counter() - t
                    t = time.perf_counter()
                self.grg.remove_mutation(mut_id, node_id)
                if instrument:
                    t_rem += time.perf_counter() - t
        self._pending_bubbles.clear()

        if not self.defer_sample_updates:
            self.flush_sample_updates()

        if instrument:
            elapsed = time.perf_counter() - t_apply_start
            self.stats["apply_bubbles_time"] += elapsed
            self.stats["get_mutation_by_id_time"] += t_get
            self.stats["add_mutation_time"] += t_add
            self.stats["remove_mutation_time"] += t_rem
            self.stats["get_mutation_by_id_calls"] += n_muts_total
            self.stats["add_mutation_calls"] += n_muts_total
            self.stats["remove_mutation_calls"] += n_muts_total
            self.stats["bubbles_created"] += n_bubbles
            self.stats["mutations_moved"] += n_muts_total

    def flush_sample_updates(self):
        if self._pending_sample_removals:
            instrument = self.instrument
            if instrument:
                t = time.perf_counter()
            current = set(self.grg.get_sample_nodes())
            current.difference_update(self._pending_sample_removals)
            self.grg.set_samples(list(current))
            self._pending_sample_removals.clear()
            if instrument:
                self.stats["flush_samples_time"] += time.perf_counter() - t

    def _clear_modified_caches(self):
        if self.instrument:
            t = time.perf_counter()
        for node_id in self._modified_nodes:
            self._mutation_cache.pop(node_id, None)
            self._pos_cache.pop(node_id, None)
            if node_id < len(self.span_cache):
                self.span_cache[node_id] = False
                self.anc_cov_cache[node_id] = False
        self._modified_nodes.clear()
        if self.instrument:
            self.stats["clear_caches_time"] += time.perf_counter() - t

    # ------------------------------------------------------------------
    # Public entry points
    # ------------------------------------------------------------------

    def _register_offspring(self, offspring_id):
        idx = self._negative_node_index.get(offspring_id)
        if idx is None:
            idx = len(self.NEGATIVE_NODE_IDS)
            self._negative_node_index[offspring_id] = idx
            self.NEGATIVE_NODE_IDS.append(offspring_id)
        return -(idx + 1)

    def recombine(self, haplotype_A, haplotype_B, breakpoint):
        instrument = self.instrument
        self._pending_bubbles.clear()

        if instrument:
            t = time.perf_counter()
            gen_v_before = self._gen_visited
        self._sync_to_grg()
        if instrument:
            self.stats["sync_to_grg_time"] += time.perf_counter() - t
            t = time.perf_counter()
        offspring_id = self.grg.make_node(negative=True)
        if instrument:
            self.stats["make_node_time"] += time.perf_counter() - t
            self.stats["make_node_calls"] += 1
        self._grow_node_arrays(self.grg.num_nodes - 1)
        self._gen_connected += 1

        if instrument:
            recurse_t = 0.0
            ts = time.perf_counter()
        self._recurse_attach(haplotype_A, offspring_id, 0, breakpoint)
        if instrument:
            recurse_t += time.perf_counter() - ts
            self.stats["recurse_attach_calls"] += 1
            self.stats["segments_processed"] += 1
            ts = time.perf_counter()
        self._recurse_attach(haplotype_B, offspring_id, breakpoint, self.genome_length)
        if instrument:
            recurse_t += time.perf_counter() - ts
            self.stats["recurse_attach_calls"] += 1
            self.stats["segments_processed"] += 1
            self.stats["recurse_attach_time"] += recurse_t

        self._apply_pending_bubbles()
        self._clear_modified_caches()

        if instrument:
            gen_v_after = self._gen_visited
            gen_set = set(range(gen_v_before + 1, gen_v_after + 1))
            v = sum(1 for x in self._visited_gen if x in gen_set)
            self.stats["visits_total"] += v
            self.stats["offspring_count"] += 1

        return self._register_offspring(offspring_id)

    def recombine_multi(self, segments):
        instrument = self.instrument
        self._pending_bubbles.clear()

        if instrument:
            t = time.perf_counter()
            gen_v_before = self._gen_visited
        self._sync_to_grg()
        if instrument:
            self.stats["sync_to_grg_time"] += time.perf_counter() - t
            t = time.perf_counter()
        offspring_id = self.grg.make_node(negative=True)
        if instrument:
            self.stats["make_node_time"] += time.perf_counter() - t
            self.stats["make_node_calls"] += 1
        self._grow_node_arrays(self.grg.num_nodes - 1)
        self._gen_connected += 1

        if instrument:
            recurse_t = 0.0
        start = 0
        for parent_id, end in segments:
            if end > start:
                if self.debug_mode:
                    print("BREAK")
                if instrument:
                    ts = time.perf_counter()
                self._recurse_attach(parent_id, offspring_id, start, end)
                if instrument:
                    recurse_t += time.perf_counter() - ts
                    self.stats["recurse_attach_calls"] += 1
                    self.stats["segments_processed"] += 1
            start = end
        if instrument:
            self.stats["recurse_attach_time"] += recurse_t

        self._apply_pending_bubbles()
        self._clear_modified_caches()

        if instrument:
            gen_v_after = self._gen_visited
            gen_set = set(range(gen_v_before + 1, gen_v_after + 1))
            v = sum(1 for x in self._visited_gen if x in gen_set)
            self.stats["visits_total"] += v
            self.stats["offspring_count"] += 1

        return self._register_offspring(offspring_id)


def simulate_grg_recombination(benchmark, recomb, bp_range, N):
    total_bp = 0
    samples = np.array(recomb.grg.get_sample_nodes())
    np.random.shuffle(samples)

    new_offspring_ids = []

    # Batch sample updates: the wholesale set_samples() at the end of this
    # function replaces samples entirely, so per-offspring sample-removal
    # work in _apply_pending_bubbles is redundant.
    prev_defer = recomb.defer_sample_updates
    recomb.defer_sample_updates = True

    try:
        for i in range(0, len(samples) - 1, 2):
            p1 = samples[i]
            p2 = samples[i + 1]

            for j in range(benchmark.num_offspring_per_couple):
                bp, num_bp = benchmark.get_breakpoints(bp_range, expected_crossovers=1.5)
                segments = recombination_intervals(p1, p2, bp, N)
                total_bp += num_bp

                offspring_id = recomb.recombine_multi(segments)
                raw_id = recomb.NEGATIVE_NODE_IDS[abs(offspring_id) - 1]
                new_offspring_ids.append(raw_id)

        new_offspring_ids.sort()
        recomb.grg.set_samples(new_offspring_ids)
        # Wholesale set_samples just made every accumulated removal moot.
        recomb._pending_sample_removals.clear()
    finally:
        recomb.defer_sample_updates = prev_defer

    return new_offspring_ids, total_bp
