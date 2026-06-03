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

    mut_counts = [len(grg.get_mutations_for_node(i, allow_sort=False)) for i in range(n)]
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
    - Up-edge walks in `_recurse_attach` skip parents already attached to
      the current offspring (`connected_gen[parent] == gen_c`). Bubbles
      created in earlier segments of the same offspring carry mutations
      confined to those earlier segments' intervals, so revisiting them
      in later segments produces only pruning_root events. Filter is
      per-offspring (gen_c bumps per recombine call), so bubbles from
      previous offspring remain reachable for legitimate reuse.
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

        # ----- Audit 1: case-counter instrumentation -----
        # Cumulative per-decision-case counts. Every visit in _recurse_attach
        # falls into exactly one of 13 matrix cells (9 standard + 4 root
        # variants); audit_check() asserts that the implementation's bubble,
        # connect, and make_node totals match what the matrix prescribes.
        # Always-on (independent of `instrument`): each increment is a single
        # dict op, dominated by surrounding work.
        self.audit = self._fresh_audit()

    @staticmethod
    def _fresh_audit():
        return {
            # row "no relevant" (Mu cap I = empty)
            'pruning':                  0,  # has_no_relevant + ancestral_disjoint
            'pruning_root':             0,  # has_no_relevant + Iu is None
            'path_compression':         0,  # has_no_relevant + full_coverage (recursed)
            'path_compression_attach':  0,  # has_no_relevant + full_coverage + fanout>1 (attached)
            'decomposition':            0,  # has_no_relevant + partial overlap
            # row "all covered" (Mu subset I)
            'direct_attach':            0,  # has_all_relevant + full_coverage
            'direct_attach_root':       0,  # has_all_relevant + Iu is None
            'bubble_strip':             0,  # has_all_relevant + ancestral_disjoint
            'bubble_split':             0,  # has_all_relevant + partial overlap
            'direct_attach_dup':        0,  # connected_gen guard skipped a duplicate edge
            # row "partial relevant"
            'bubble_fill':              0,  # has_partial_relevant + full_coverage
            'bubble_strip_partial':     0,  # has_partial_relevant + ancestral_disjoint
            'bubble_split_partial':     0,  # has_partial_relevant + partial overlap
            'bubble_strip_partial_rt':  0,  # has_partial_relevant + Iu is None
            # informational skips (not in matrix)
            'skip_empty_interval':      0,  # L >= R at top of loop
            'skip_already_visited':     0,  # gen_v guard fired
            'skip_empty_trim':          0,  # newL >= newR after trim
            # primitives that should reconcile against case sums
            'visits':                   0,  # entries past both skip guards
            'extract_bubble_calls':     0,
            'connect_calls_in_attach':  0,  # connects fired inside _recurse_attach
            'connect_calls_in_extract': 0,  # connects fired inside _extract_bubble (2/call)
            'make_node_calls':          0,  # offspring + bubble nodes created
            'recombine_calls':          0,  # recombine() + recombine_multi() entries
            'recurse_attach_calls':     0,
        }

    @staticmethod
    def _fresh_stats():
        return {
            # Phase-level wallclock totals (seconds)
            "recurse_attach_time": 0.0,
            "apply_bubbles_time": 0.0,
            "sync_to_grg_time": 0.0,
            "clear_caches_time": 0.0,
            "flush_samples_time": 0.0,
            "sort_mutations_time": 0.0,
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
        self.audit = self._fresh_audit()

    def audit_reset(self):
        """Zero only the audit counters (leaves stats untouched)."""
        self.audit = self._fresh_audit()

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
            mut_ids = self.grg.get_mutations_for_node(node_id, allow_sort=False)
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
        audit = self.audit
        if instrument:
            t = time.perf_counter()
        bubble_id = self.grg.make_node()
        if instrument:
            self.stats["make_node_time"] += time.perf_counter() - t
            self.stats["make_node_calls"] += 1
        audit['make_node_calls'] += 1

        self._grow_node_arrays(bubble_id)

        if instrument:
            t = time.perf_counter()
        self.grg.connect(bubble_id, node_id)
        self.grg.connect(bubble_id, -offspring_id)
        if instrument:
            self.stats["connect_time"] += time.perf_counter() - t
            self.stats["connect_calls"] += 2
        audit['connect_calls_in_extract'] += 2
        audit['extract_bubble_calls'] += 1

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
        audit = self.audit  # hoisted for inner-loop speed

        audit['recurse_attach_calls'] += 1

        neg_offspring = -offspring_id
        stack = [(root_id, L0, R0)]
        stack_append = stack.append
        stack_pop = stack.pop

        while stack:
            node_id, L, R = stack_pop()

            if L >= R:
                audit['skip_empty_interval'] += 1
                continue
            if visited_gen[node_id] == gen_v:
                audit['skip_already_visited'] += 1
                continue
            visited_gen[node_id] = gen_v
            audit['visits'] += 1

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
                # Root: only own muts matter. The 4 sub-cases are mutually
                # exclusive; structured as if/elif/else so audit tagging is
                # one-per-visit.
                if has_all_relevant:
                    if connected_gen[node_id] != gen_c:
                        grg_connect(node_id, neg_offspring)
                        connected_gen[node_id] = gen_c
                        pending_sample_removals_add(node_id)
                        audit['direct_attach_root'] += 1
                        audit['connect_calls_in_attach'] += 1
                    else:
                        audit['direct_attach_dup'] += 1
                elif has_partial_relevant:
                    rel_mut_ids = [m[0] for m in mutation_cache[node_id][left:right]]
                    bubble_id = extract_bubble(node_id, rel_mut_ids, offspring_id, (L, R))
                    connected_gen[bubble_id] = gen_c
                    audit['bubble_strip_partial_rt'] += 1
                else:
                    audit['pruning_root'] += 1
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
                    audit['direct_attach'] += 1
                    audit['connect_calls_in_attach'] += 1
                else:
                    audit['direct_attach_dup'] += 1
                continue

            if has_no_relevant and ancestral_disjoint:
                audit['pruning'] += 1
                continue

            if has_no_relevant:
                parents = get_up_edges_cached(node_id)
                # Early-attach optimization: if a multi-parent empty
                # intermediary fully covers [L, R), attach the offspring
                # to it directly instead of walking up to its ancestors.
                # Saves K-1 edges per chain (K = number of attaches the
                # avoided recursion would have made). Gated on:
                #   num_all == 0      -- has_no_relevant also matches nodes
                #                        whose mutations all sit outside
                #                        [L, R); attaching offspring there
                #                        would transitively pull in those
                #                        out-of-segment muts (correctness).
                #   full_coverage     -- Iu(node_id) ⊂ [L, R), so the
                #                        offspring inherits exactly the
                #                        right mutations via this path.
                #   len(parents) > 1  -- at fanout-1 the chain compresses
                #                        to one edge anyway; attaching here
                #                        is strictly worse (same edges,
                #                        +1 hop forever).
                #
                # Multitree preservation: no extra work is required, by the
                # same argument that lets direct_attach skip recursion. In
                # a multitree input, every ancestor of node_id reaches the
                # haplotype only via node_id, so sibling paths in this
                # `_recurse_attach` call cannot reach those ancestors and
                # cannot fire any attach branch on them. Bubbles created
                # mid-call are roots (no parents), so they don't reintroduce
                # multi-route paths. The frontier remains an antichain.
                if num_all == 0 and full_coverage and len(parents) > 1:
                    if connected_gen[node_id] != gen_c:
                        grg_connect(node_id, neg_offspring)
                        connected_gen[node_id] = gen_c
                        pending_sample_removals_add(node_id)
                        audit['path_compression_attach'] += 1
                        audit['connect_calls_in_attach'] += 1
                    else:
                        audit['direct_attach_dup'] += 1
                    continue
                # Inline max/min — Python's builtins are surprisingly slow
                # at this call rate (showed up at ~50ms in the profile).
                newL = L if L > Iu0 else Iu0
                newR = R if R < Iu1 else Iu1
                if full_coverage:
                    audit['path_compression'] += 1
                else:
                    audit['decomposition'] += 1
                if newL >= newR:
                    audit['skip_empty_trim'] += 1
                    continue
                # Skip parents already attached to this offspring (typically
                # bubbles created in earlier segments). Their mutations live
                # in disjoint intervals so re-descending them only produces
                # pruning_root events.
                for parent in reversed(parents):
                    if connected_gen[parent] != gen_c:
                        stack_append((parent, newL, newR))
                continue

            # Partial / all relevant, not full coverage -> bubble + maybe recurse.
            rel_mut_ids = [m[0] for m in mutation_cache[node_id][left:right]]
            bubble_id = extract_bubble(node_id, rel_mut_ids, offspring_id, (L, R))
            connected_gen[bubble_id] = gen_c

            # Classify exactly one of the 5 non-root bubble cells.
            if has_all_relevant:
                if ancestral_disjoint:
                    audit['bubble_strip'] += 1
                else:
                    audit['bubble_split'] += 1
            else:  # has_partial_relevant
                if full_coverage:
                    audit['bubble_fill'] += 1
                elif ancestral_disjoint:
                    audit['bubble_strip_partial'] += 1
                else:
                    audit['bubble_split_partial'] += 1

            if not ancestral_disjoint:
                newL = L if L > Iu0 else Iu0
                newR = R if R < Iu1 else Iu1
                if newL >= newR:
                    audit['skip_empty_trim'] += 1
                    continue
                # Same offspring-local skip: already-attached parents have
                # nothing useful to contribute to subsequent segments.
                for parent in reversed(get_up_edges_cached(node_id)):
                    if connected_gen[parent] != gen_c:
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
        self.audit['recombine_calls'] += 1
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
        self.audit['make_node_calls'] += 1
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
        self.audit['recombine_calls'] += 1
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
        self.audit['make_node_calls'] += 1
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

    # ------------------------------------------------------------------
    # Audit 1: invariant checks and reporting
    # ------------------------------------------------------------------

    def audit_check(self, audit=None, raise_on_fail=True):
        """
        Verify implementation matches algorithm spec on cumulative counts.

        Pass an `audit` dict to check an aggregated/saved snapshot;
        otherwise checks `self.audit`. Raises AssertionError on first
        failure if `raise_on_fail`. Returns dict of {name: result}.

        Invariants:
          (1) extract_bubble_calls == sum of 6 bubble-case counters.
              Failure => bubble was created outside the matrix, or a
              matrix branch ran without calling _extract_bubble.
          (2) total connect calls == 2*extract_bubble_calls + firing
              direct attaches. Failure => connect() was invoked outside
              the two prescribed paths.
          (3) make_node_calls == recombine_calls + extract_bubble_calls.
              Failure => a node was created outside offspring/bubble paths.
        """
        a = audit if audit is not None else self.audit

        bubble_cases = (a['bubble_strip'] + a['bubble_split'] +
                        a['bubble_fill'] + a['bubble_strip_partial'] +
                        a['bubble_split_partial'] + a['bubble_strip_partial_rt'])
        # direct_attach_dup did NOT call connect (gen_c guard fired),
        # so it is excluded from the connect identity. path_compression_attach
        # DOES call connect (one per firing), so it joins the direct-firing
        # term that balances against connect_calls_in_attach.
        direct_firing = (a['direct_attach'] + a['direct_attach_root']
                         + a['path_compression_attach'])
        total_connects = a['connect_calls_in_attach'] + a['connect_calls_in_extract']
        expected_connects = 2 * a['extract_bubble_calls'] + direct_firing
        expected_make_nodes = a['recombine_calls'] + a['extract_bubble_calls']

        results = {
            'bubble_identity': {
                'lhs': a['extract_bubble_calls'], 'rhs': bubble_cases,
                'pass': a['extract_bubble_calls'] == bubble_cases,
                'desc': 'extract_bubble_calls == sum of 6 bubble-case counters',
            },
            'connect_identity': {
                'lhs': total_connects, 'rhs': expected_connects,
                'pass': total_connects == expected_connects,
                'desc': 'total connect calls == 2*bubbles + firing direct attaches',
            },
            'make_node_identity': {
                'lhs': a['make_node_calls'], 'rhs': expected_make_nodes,
                'pass': a['make_node_calls'] == expected_make_nodes,
                'desc': 'make_node_calls == recombine_calls + bubbles',
            },
        }

        if raise_on_fail:
            for name, r in results.items():
                if not r['pass']:
                    raise AssertionError(
                        f"audit_check FAIL [{name}]: {r['desc']} -- "
                        f"lhs={r['lhs']}, rhs={r['rhs']}, "
                        f"delta={r['lhs'] - r['rhs']}"
                    )
        return results

    def audit_summary(self, audit=None, header=None):
        """Pretty-print the case histogram + identity check.

        Pass `audit` to summarise an aggregated/saved snapshot; otherwise
        summarises `self.audit`. `header` is an optional title prefix.
        """
        a = audit if audit is not None else self.audit
        bubble_cases = (a['bubble_strip'] + a['bubble_split'] +
                        a['bubble_fill'] + a['bubble_strip_partial'] +
                        a['bubble_split_partial'] + a['bubble_strip_partial_rt'])
        direct_firing = (a['direct_attach'] + a['direct_attach_root']
                         + a['path_compression_attach'])
        total_decisions = (
            a['pruning'] + a['pruning_root']
            + a['path_compression'] + a['path_compression_attach']
            + a['decomposition']
            + a['direct_attach'] + a['direct_attach_root'] + a['direct_attach_dup']
            + bubble_cases
        )

        def pct(n):
            return (100.0 * n / total_decisions) if total_decisions else 0.0

        lines = []
        title = "AUDIT 1 -- decision case histogram"
        if header:
            title = f"{title}  ({header})"
        lines.append("=" * 64)
        lines.append(title)
        lines.append("=" * 64)
        lines.append(f"{'case':<30s} {'count':>14s}  {'%':>6s}")
        lines.append("-" * 64)
        for k in ('pruning', 'pruning_root',
                  'path_compression', 'path_compression_attach',
                  'decomposition',
                  'direct_attach', 'direct_attach_root', 'direct_attach_dup',
                  'bubble_strip', 'bubble_split',
                  'bubble_fill',
                  'bubble_strip_partial', 'bubble_split_partial',
                  'bubble_strip_partial_rt'):
            lines.append(f"{k:<30s} {a[k]:>14d}  {pct(a[k]):>5.1f}%")
        lines.append("-" * 64)
        lines.append(f"{'TOTAL DECISIONS':<30s} {total_decisions:>14d}")
        lines.append("")
        lines.append(f"{'visits (post-skip-guards)':<30s} {a['visits']:>14d}")
        lines.append(f"{'skip_empty_interval':<30s} {a['skip_empty_interval']:>14d}")
        lines.append(f"{'skip_already_visited':<30s} {a['skip_already_visited']:>14d}")
        lines.append(f"{'skip_empty_trim':<30s} {a['skip_empty_trim']:>14d}")
        lines.append("")
        lines.append("Primitives:")
        lines.append(f"{'recombine_calls':<30s} {a['recombine_calls']:>14d}")
        lines.append(f"{'recurse_attach_calls':<30s} {a['recurse_attach_calls']:>14d}")
        lines.append(f"{'extract_bubble_calls':<30s} {a['extract_bubble_calls']:>14d}")
        lines.append(f"{'make_node_calls':<30s} {a['make_node_calls']:>14d}")
        lines.append(f"{'connect (in _recurse_attach)':<30s} {a['connect_calls_in_attach']:>14d}")
        lines.append(f"{'connect (in _extract_bubble)':<30s} {a['connect_calls_in_extract']:>14d}")
        lines.append("")
        if a['recombine_calls']:
            lines.append("Per-recombine averages:")
            lines.append(f"  bubbles / offspring     {a['extract_bubble_calls']/a['recombine_calls']:>10.3f}")
            lines.append(f"  visits / offspring      {a['visits']/a['recombine_calls']:>10.3f}")
            lines.append(f"  segments / offspring    {a['recurse_attach_calls']/a['recombine_calls']:>10.3f}")
            lines.append(f"  edges added / offspring {(2*bubble_cases + direct_firing)/a['recombine_calls']:>10.3f}")
        lines.append("")
        lines.append("Identities (via audit_check):")
        results = self.audit_check(audit=a, raise_on_fail=False)
        for r in results.values():
            mark = "OK" if r['pass'] else "FAIL"
            lines.append(f"  [{mark}] {r['desc']}: lhs={r['lhs']} rhs={r['rhs']}")
        lines.append("=" * 64)
        print("\n".join(lines))


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

        # Compact soft-deleted mutation entries left behind by remove_mutation
        # (now O(1) with the new API). Without this, get_mutations_for_node
        # walks past the dead entries on every read at heavily-modified nodes.
        if recomb.instrument:
            t_sort = time.perf_counter()
        recomb.grg.sort_mutations()
        if recomb.instrument:
            recomb.stats["sort_mutations_time"] += time.perf_counter() - t_sort
        # sort_mutations renumbers MutationIds (by (position, allele)), so any
        # cached (mut_id, position) pairs from this generation are now stale.
        # Position-derived caches (span_cache, anc_cov_cache, _up_edges_cache)
        # are unaffected and survive. Lazy repopulation on first touch in the
        # next generation.
        recomb._mutation_cache.clear()
        recomb._pos_cache.clear()
    finally:
        recomb.defer_sample_updates = prev_defer

    return new_offspring_ids, total_bp
