"""
grg_recombination.py

Optimized core module containing the Non-duplication GRG recombination algorithm.
"""

import bisect
import numpy as np

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
    """

    debug_mode = False

    def __init__(self, grg):
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
        bubble_id = self.grg.make_node()
        self._grow_node_arrays(bubble_id)

        self.grg.connect(bubble_id, node_id)
        self.grg.connect(bubble_id, -offspring_id)

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
        for bubble_op in self._pending_bubbles:
            node_id = bubble_op['node_id']
            bubble_id = bubble_op['bubble_id']
            for mut_id in bubble_op['relevant_mut_ids']:
                mut = self.grg.get_mutation_by_id(mut_id)
                self.grg.add_mutation(mut, bubble_id)
                self.grg.remove_mutation(mut_id, node_id)
        self._pending_bubbles.clear()

        if not self.defer_sample_updates:
            self.flush_sample_updates()

    def flush_sample_updates(self):
        if self._pending_sample_removals:
            current = set(self.grg.get_sample_nodes())
            current.difference_update(self._pending_sample_removals)
            self.grg.set_samples(list(current))
            self._pending_sample_removals.clear()

    def _clear_modified_caches(self):
        for node_id in self._modified_nodes:
            self._mutation_cache.pop(node_id, None)
            self._pos_cache.pop(node_id, None)
            if node_id < len(self.span_cache):
                self.span_cache[node_id] = False
                self.anc_cov_cache[node_id] = False
        self._modified_nodes.clear()

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
        self._pending_bubbles.clear()
        self._sync_to_grg()
        offspring_id = self.grg.make_node(negative=True)
        self._grow_node_arrays(self.grg.num_nodes - 1)
        self._gen_connected += 1

        self._recurse_attach(haplotype_A, offspring_id, 0, breakpoint)
        self._recurse_attach(haplotype_B, offspring_id, breakpoint, self.genome_length)

        self._apply_pending_bubbles()
        self._clear_modified_caches()

        return self._register_offspring(offspring_id)

    def recombine_multi(self, segments):
        self._pending_bubbles.clear()
        self._sync_to_grg()
        offspring_id = self.grg.make_node(negative=True)
        self._grow_node_arrays(self.grg.num_nodes - 1)
        self._gen_connected += 1

        start = 0
        for parent_id, end in segments:
            if end > start:
                if self.debug_mode:
                    print("BREAK")
                self._recurse_attach(parent_id, offspring_id, start, end)
            start = end

        self._apply_pending_bubbles()
        self._clear_modified_caches()

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
