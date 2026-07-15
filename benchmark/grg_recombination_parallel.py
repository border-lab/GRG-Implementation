"""
grg_recombination_parallel.py

Node-Aggregated Parallel Non-Duplication Recombination, as specified in the
"Parallel Non-Duplication Recombination" writeup (sections 1.1-1.11). This is
the parallel counterpart to `grg_recombination.NonDuplicationRecombination`
(the sequential "Insertion method" / bubble-extraction algorithm) -- same
non-duplication invariant, same offspring semantics, different traversal
strategy.

-------------------------------------------------------------------------
Two phases (writeup section 1)
-------------------------------------------------------------------------
1. Parallel discovery (`Discover`, section 1.4). Read-only w.r.t. the old
   graph. For each inherited segment (h, I), walks the ancestral cone of h
   upward through real parent edges, using the ancestor-span cache for
   pruning/whole-cone-attach decisions, and otherwise records explicit
   per-node demand. Threads only ever *submit* further work and return --
   they never block waiting on a child -- so nested fan-out cannot deadlock
   the executor.
2. Deterministic commit (`_commit`, section 1.5-1.6). The only phase that
   mutates the graph. Groups demand by old node, computes the actually
   uncovered mutation set per node (subtracting whatever direct-attach
   already transmits), and applies at most one bubble split per old node
   for this offspring plus the canonicalized direct-attach edges.

-------------------------------------------------------------------------
Mapping from writeup notation to this implementation
-------------------------------------------------------------------------
- Phi(u) (section 1.1, "ancestor span") is exactly the base class's
  `span_cache[u]` (own mutations union parents' span_cache), read via the
  inherited `_get_node_and_ancestor_span`. `span_cache[u]` is stored as an
  inclusive (min_pos, max_pos) pair; Phi(u) as a half-open interval is
  `[span[0], span[1] + 1)`.
- Edges in pygrgl carry no interval label (`grg.connect(parent, child)`
  takes no interval argument), so `Jp,u` (writeup 1.1) is always the full
  genome for every edge -- `Ip = J \u2229 Jp,u` in Algorithm 2 line 18
  degenerates to `Ip = J`. Pruning on the way up still happens via each
  parent's own Phi(p) check (Algorithm 2 line 19), so there is no loss of
  precision, only a slightly different traversal footprint than a graph
  with genuinely interval-labeled edges would have.
- `Attach(u, w, I)` (writeup 1.2) is realized as a plain `connect(u, -w)`;
  since pygrgl edges are unlabeled, "carrying interval I" is conceptual --
  the interval is implied by which mutations live in u's ancestral cone,
  which pruning has already confirmed lies inside I.
- For a single offspring w, `CovA(u, m)` (writeup 1.5) is either `{w}` or
  `\u2205`: `{w}` exactly when u ended up in the canonicalized attach set A
  (the offspring already receives the *entire* cone of u, including every
  mutation stored at u). So `Su = Mu \u2229 R[u]` when u is *not* attached,
  and `Su = \u2205` when u *is* attached (section 1.7's reasoning, applied
  to the direct-attach case rather than duplicate bubble requests). This is
  why nodes in the attach set are skipped entirely during commit.

-------------------------------------------------------------------------
Concurrency model and its limits
-------------------------------------------------------------------------
- Discovery genuinely runs on a `ThreadPoolExecutor`. Per-node `Seen[u]`/
  `R[u]` read-modify-write updates (Algorithm 2 lines 12-16) are guarded by
  a small pool of striped `threading.Lock`s (one lock per `node_id % S`
  bucket) rather than one lock per node, trading a little contention for
  O(1) memory instead of O(V) lock objects.
- The cache-population fallbacks inherited from `NonDuplicationRecombination`
  (`_get_node_and_ancestor_span`, `_get_up_edges_cached`,
  `_get_node_mutations`) are *not* additionally locked here. Two threads
  racing on the same cold cache entry can both redundantly recompute it,
  but since the computation is a pure function of already-committed graph
  state, both writes are identical and idempotent -- a benign race, not a
  correctness hazard. Introducing per-entry locks would only slow down the
  common (already-warm, from `_build_ancestral_caches`) path.
- Commit deliberately does *not* fan the actual graph-mutating pygrgl calls
  (`make_node`, `connect`, `add_mutation`, `remove_mutation`) out across
  threads: pygrgl's mutable-GRG object is a single shared native structure
  with no documented thread-safety guarantee for concurrent writers. The
  per-node *independence* the writeup describes for commit (section 1.5,
  "each old node is assigned to exactly one worker") is instead exploited
  by parallelizing only the pure-Python `Su` computation (bisecting each
  node's cached mutation positions against its recorded demand) on the
  executor, then applying every mutation/edge change from a single thread
  in a deterministic (sorted-by-node-id) order.
- Because CPython holds the GIL while executing Python bytecode, real
  wall-clock speedup from the discovery phase depends on whether pygrgl's
  C++ read accessors (`get_up_edges`, `get_mutations_for_node`, ...)
  release the GIL during the call. This module is written to be *correct*
  under genuine concurrent execution regardless; whether it is also
  *faster* than the sequential implementation is an empirical question for
  a given pygrgl build.

-------------------------------------------------------------------------
Consequence: at most one bubble per old node per offspring
-------------------------------------------------------------------------
Because `R[u]` aggregates demand from *every* segment of the current
offspring before commit runs (writeup section 1, "for each old node u, it
builds one interval-set Ru"), a single node contributes at most one bubble
per offspring here -- unlike `NonDuplicationRecombination._recurse_attach`,
which resets its `visited` generation once per *segment* and so can create
several bubbles at the same node across a multi-segment offspring. This
also means the sequential algorithm's 13-cell decision matrix collapses to
three discovery outcomes (prune / whole-cone attach / record demand) plus
one commit-side split -- the extra cases in the sequential matrix are all
about *when to stop recursing upward*, not about *what mutations move*.

The batched multi-offspring extension (writeup section 1.8, consumer
signatures) is not implemented here -- `simulate_grg_recombination` (see
`grg_recombination.py`) still processes one offspring per `recombine_multi`
call, matching the sequential algorithm's usage pattern. Only the
single-offspring algorithm (sections 1.1-1.7, 1.9-1.11) is realized.
"""

import bisect
import concurrent.futures
import os
import threading

from grg_recombination import NonDuplicationRecombination


# ---------------------------------------------------------------------------
# Canonical interval-set algebra
#
# All interval-sets in this module are lists of half-open (start, end)
# tuples, sorted by start and pairwise disjoint ("canonical" per writeup
# section 1.1). Sizes are small in practice (a handful of breakpoints per
# offspring), so these favor simplicity over asymptotic optimality.
# ---------------------------------------------------------------------------

def _normalize_intervals(intervals):
    """Sort + merge an arbitrary list of (start, end) tuples into canonical
    disjoint form. Drops empty/invalid (start >= end) tuples."""
    ivs = sorted(iv for iv in intervals if iv[0] < iv[1])
    if not ivs:
        return []
    merged = [[ivs[0][0], ivs[0][1]]]
    for s, e in ivs[1:]:
        top = merged[-1]
        if s <= top[1]:
            if e > top[1]:
                top[1] = e
        else:
            merged.append([s, e])
    return [tuple(iv) for iv in merged]


def _interval_union(a, b):
    """Canonical union of two (already-canonical) interval-sets."""
    if not a:
        return list(b)
    if not b:
        return list(a)
    return _normalize_intervals(list(a) + list(b))


def _interval_diff(a, b):
    """a \\ b, for canonical sorted-disjoint interval-sets a, b."""
    if not a or not b:
        return list(a)
    result = []
    for s, e in a:
        cur = s
        for bs, be in b:
            if be <= cur:
                continue
            if bs >= e:
                break
            if bs > cur:
                result.append((cur, bs))
            if be > cur:
                cur = be
            if cur >= e:
                break
        if cur < e:
            result.append((cur, e))
    return result


def _interval_overlaps_range(intervals, lo, hi):
    """True iff the interval-set `intervals` overlaps the half-open range
    [lo, hi)."""
    for s, e in intervals:
        if e > lo and s < hi:
            return True
    return False


def _interval_covers_range(intervals, lo, hi):
    """True iff the union of (canonical, sorted) `intervals` fully covers
    the half-open range [lo, hi), i.e. no gaps."""
    if lo >= hi:
        return True
    cur = lo
    for s, e in intervals:
        if e <= cur:
            continue
        if s > cur:
            return False
        if e > cur:
            cur = e
        if cur >= hi:
            return True
    return cur >= hi


# ---------------------------------------------------------------------------
# Recursive fan-out helper
# ---------------------------------------------------------------------------

class _FanOut:
    """Tracks outstanding recursively-submitted tasks on a shared executor.

    Every worker only ever submits child tasks and returns -- it never
    blocks waiting on a child's result -- so nested `submit()` calls cannot
    starve/deadlock a bounded thread pool. `wait()` blocks the *calling*
    (non-pool) thread until every transitively-submitted task has finished,
    then re-raises the first exception observed, if any.

    The counter starts at 1 (a "phantom" reference held by the initiating
    thread) rather than 0. Without it, a worker could race ahead and finish
    the *first* top-level task -- dropping the real count back to 0 -- while
    the initiating thread is still partway through submitting the second,
    third, ... top-level task, firing `_done` and letting `wait()` return
    before all the real work even exists. `close_initial()` releases the
    phantom reference once every top-level `submit()` call has been issued;
    only then can the count legitimately reach 0.
    """

    __slots__ = ("_executor", "_lock", "_count", "_done", "_exc")

    def __init__(self, executor):
        self._executor = executor
        self._lock = threading.Lock()
        self._count = 1
        self._done = threading.Event()
        self._exc = None

    def submit(self, fn, *args):
        with self._lock:
            self._count += 1
        self._executor.submit(self._run, fn, args)

    def _run(self, fn, args):
        try:
            fn(*args)
        except BaseException as exc:  # noqa: BLE001 - propagate to waiter
            with self._lock:
                if self._exc is None:
                    self._exc = exc
        finally:
            self._release()

    def _release(self):
        with self._lock:
            self._count -= 1
            count = self._count
        if count == 0:
            self._done.set()

    def close_initial(self):
        """Release the phantom initial reference. Call exactly once, after
        every top-level `submit()` has been issued, before `wait()`."""
        self._release()

    def wait(self):
        self._done.wait()
        if self._exc is not None:
            raise self._exc


class ParallelNonDuplicationRecombination(NonDuplicationRecombination):
    """
    Node-Aggregated Parallel Non-Duplication Recombination (single-offspring
    case). See module docstring for the full correspondence to the writeup
    and the concurrency trade-offs.

    Public interface intentionally mirrors `NonDuplicationRecombination`
    (`recombine`, `recombine_multi`, `audit_check`, ...) so it is a drop-in
    replacement everywhere `simulate_grg_recombination` is used.
    """

    def __init__(self, grg, instrument=False, max_workers=None):
        super().__init__(grg, instrument=instrument)

        if max_workers is None:
            max_workers = os.cpu_count() or 4
        self.max_workers = max_workers
        self._executor = concurrent.futures.ThreadPoolExecutor(
            max_workers=max_workers, thread_name_prefix="grg-discover"
        )

        # Lock striping for Seen[u]/R[u] read-modify-write (Algorithm 2,
        # lines 12-16). A fixed-size pool avoids O(V) lock objects on large
        # graphs; correctness only requires that concurrent claims on the
        # *same* node_id serialize, which striping still guarantees.
        self._num_lock_stripes = 256
        self._node_locks = [threading.Lock() for _ in range(self._num_lock_stripes)]

        # Guards the discovery-only diagnostic counters below (bumped from
        # worker threads); the commit-phase counters are main-thread-only
        # and need no lock.
        self._audit_lock = threading.Lock()

    def _lock_for(self, node_id):
        return self._node_locks[node_id % self._num_lock_stripes]

    @staticmethod
    def _fresh_audit():
        # Extend the base (sequential) audit dict rather than replacing it:
        # `_extract_bubble` (inherited, used unmodified in `_commit`) bumps
        # 'make_node_calls', 'connect_calls_in_extract', and
        # 'extract_bubble_calls' by the base class's names, so those keys
        # must exist. The sequential decision-matrix keys (bubble_strip,
        # direct_attach, ...) are simply never touched here.
        audit = NonDuplicationRecombination._fresh_audit()
        audit.update({
            # Discovery-phase counters (bumped under _audit_lock, only when
            # self.instrument is True -- see _bump).
            'discover_pruned':          0,  # Phi(u) ∩ I == ∅
            'discover_attach_whole':    0,  # Phi(u) ⊆ I -> whole-cone attach
            'discover_claim_empty':     0,  # J = I \ Seen[u] == ∅
            'discover_demand_recorded': 0,  # J ≠ ∅ -> recorded into R[u]
            # Commit-phase counters (main-thread only, always on).
            'nodes_with_demand_total':  0,  # sum of |{u : R[u] != ∅}| per offspring
            'nodes_attached_total':     0,  # sum of |canonicalized A| per offspring
            'nodes_covered_by_attach':  0,  # u in both Usplit and A -> Su forced ∅
            'mutations_moved':          0,
        })
        return audit

    def _bump(self, key):
        """Discovery-phase counter increment; gated on `instrument` since
        it is the only place multiple threads touch `self.audit`
        concurrently and taking a lock on every visit isn't free."""
        if not self.instrument:
            return
        with self._audit_lock:
            self.audit[key] += 1

    def shutdown(self, wait=True):
        self._executor.shutdown(wait=wait)

    def __del__(self):
        try:
            self._executor.shutdown(wait=False)
        except Exception:
            pass

    # ------------------------------------------------------------------
    # Phase 1: parallel discovery (writeup Algorithm 2)
    # ------------------------------------------------------------------

    def _discover_task(self, node_id, intervals, state):
        """Discover(u=node_id, w=state['w'], I=intervals, A=state['attach'],
        R=state['demand'], Seen=state['seen']).

        Read-only w.r.t. the old graph: only ever writes into `state`'s
        thread-shared structures (attach/demand/seen), guarded as described
        in the module docstring.
        """
        if not intervals:
            return

        span = self._get_node_and_ancestor_span(node_id)  # Phi(u), inclusive
        if span is None:
            self._bump('discover_pruned')
            return

        lo = span[0]
        hi = span[1] + 1  # half-open

        if not _interval_overlaps_range(intervals, lo, hi):
            self._bump('discover_pruned')
            return

        if _interval_covers_range(intervals, lo, hi):
            with state['attach_lock']:
                state['attach'].add(node_id)
            self._bump('discover_attach_whole')
            return

        lock = self._lock_for(node_id)
        with lock:
            prev_seen = state['seen'].get(node_id, ())
            claimed = _interval_diff(intervals, prev_seen)
            if not claimed:
                self._bump('discover_claim_empty')
                return
            state['seen'][node_id] = _interval_union(prev_seen, claimed)
            prev_demand = state['demand'].get(node_id, ())
            state['demand'][node_id] = _interval_union(prev_demand, claimed)

        self._bump('discover_demand_recorded')

        fanout = state['fanout']
        for parent in self._get_up_edges_cached(node_id):
            p_span = self._get_node_and_ancestor_span(parent)
            if p_span is None:
                continue
            # Ip = J ∩ Jp,u = claimed (pygrgl edges carry no interval label);
            # pre-check Phi(p) here purely to avoid scheduling dead work --
            # _discover_task would reach the same conclusion itself.
            if _interval_overlaps_range(claimed, p_span[0], p_span[1] + 1):
                fanout.submit(self._discover_task, parent, claimed, state)

    def _discover_all(self, segments, offspring_id):
        state = {
            'seen': {},
            'demand': {},
            'attach': set(),
            'attach_lock': threading.Lock(),
            'w': offspring_id,
        }
        fanout = _FanOut(self._executor)
        state['fanout'] = fanout

        start = 0
        for parent_id, end in segments:
            if end > start:
                fanout.submit(self._discover_task, parent_id, [(start, end)], state)
            start = end

        fanout.close_initial()
        fanout.wait()
        return state

    # ------------------------------------------------------------------
    # Phase 2: deterministic commit (writeup Algorithm 3/4)
    # ------------------------------------------------------------------

    def _mutation_ids_in_intervals(self, node_id, intervals):
        """Mu ∩ (union of intervals), as a flat list of mutation ids."""
        self._get_node_mutations(node_id)
        positions = self._pos_cache[node_id]
        if not positions:
            return []
        mutation_cache = self._mutation_cache[node_id]
        result = []
        for lo, hi in intervals:
            left = bisect.bisect_left(positions, lo)
            right = bisect.bisect_left(positions, hi)
            if right > left:
                result.extend(m[0] for m in mutation_cache[left:right])
        return result

    def _commit(self, offspring_id, state):
        neg_offspring = -offspring_id
        attach_set = state['attach']
        demand = state['demand']
        audit = self.audit

        audit['nodes_attached_total'] += len(attach_set)
        for node_id in sorted(attach_set):
            self.grg.connect(node_id, neg_offspring)
            self._pending_sample_removals.add(node_id)
            audit['connect_calls_in_attach'] += 1

        # Su = {m in Mu ∩ R[u] : w not in CovA(u, m)}. For a single
        # offspring, CovA(u, m) = {w} for every m in Mu when u is in the
        # canonicalized attach set (the whole cone -- including all of Mu --
        # already flows to w directly), else ∅. So nodes in attach_set
        # contribute nothing to Su regardless of what else was recorded in
        # R[u] (writeup section 1.7's covered-request argument, applied to
        # direct attach instead of a duplicate bubble request).
        usplit = sorted(u for u, ivs in demand.items() if ivs and u not in attach_set)
        audit['nodes_with_demand_total'] += len(demand)
        audit['nodes_covered_by_attach'] += sum(1 for u in demand if u in attach_set)

        if usplit:
            su_lists = list(self._executor.map(
                self._mutation_ids_in_intervals, usplit, [demand[u] for u in usplit]
            ))
        else:
            su_lists = []

        for node_id, rel_mut_ids in zip(usplit, su_lists):
            if not rel_mut_ids:
                continue
            ivs = demand[node_id]
            bubble_interval = (ivs[0][0], ivs[-1][1])
            self._extract_bubble(node_id, rel_mut_ids, offspring_id, bubble_interval)
            audit['mutations_moved'] += len(rel_mut_ids)

    # ------------------------------------------------------------------
    # Public entry points (mirror NonDuplicationRecombination's interface)
    # ------------------------------------------------------------------

    def recombine_multi(self, segments):
        self.audit['recombine_calls'] += 1
        self._pending_bubbles.clear()

        self._sync_to_grg()
        offspring_id = self.grg.make_node(negative=True)
        self.audit['make_node_calls'] += 1
        self._grow_node_arrays(self.grg.num_nodes - 1)

        state = self._discover_all(segments, offspring_id)
        self._commit(offspring_id, state)

        self._apply_pending_bubbles()
        self._clear_modified_caches()

        return self._register_offspring(offspring_id)

    def recombine(self, haplotype_A, haplotype_B, breakpoint):
        return self.recombine_multi(
            [(haplotype_A, breakpoint), (haplotype_B, self.genome_length)]
        )

    # ------------------------------------------------------------------
    # Audit 1 (parallel variant): invariant checks and reporting
    # ------------------------------------------------------------------

    def audit_check(self, audit=None, raise_on_fail=True):
        """
        Verify implementation matches algorithm spec on cumulative counts.

        Invariants:
          (1) make_node_calls == recombine_calls + extract_bubble_calls.
              Failure => a node was created outside offspring/bubble paths.
          (2) connect_calls_in_extract == 2 * extract_bubble_calls.
              Failure => _extract_bubble's two-connect contract was broken.
          (3) connect_calls_in_attach == nodes_attached_total.
              Failure => connect() fired outside the canonicalized attach
              loop, or fired more/fewer times than canonicalized attaches.
          (4) total connects == nodes_attached_total + 2*extract_bubble_calls.
              Failure => connect() was invoked outside the two prescribed
              commit-phase paths.
        """
        a = audit if audit is not None else self.audit

        expected_make_nodes = a['recombine_calls'] + a['extract_bubble_calls']
        expected_extract_connects = 2 * a['extract_bubble_calls']
        total_connects = a['connect_calls_in_attach'] + a['connect_calls_in_extract']
        expected_total_connects = a['nodes_attached_total'] + expected_extract_connects

        results = {
            'make_node_identity': {
                'lhs': a['make_node_calls'], 'rhs': expected_make_nodes,
                'pass': a['make_node_calls'] == expected_make_nodes,
                'desc': 'make_node_calls == recombine_calls + extract_bubble_calls',
            },
            'bubble_connect_identity': {
                'lhs': a['connect_calls_in_extract'], 'rhs': expected_extract_connects,
                'pass': a['connect_calls_in_extract'] == expected_extract_connects,
                'desc': 'connect_calls_in_extract == 2 * extract_bubble_calls',
            },
            'attach_connect_identity': {
                'lhs': a['connect_calls_in_attach'], 'rhs': a['nodes_attached_total'],
                'pass': a['connect_calls_in_attach'] == a['nodes_attached_total'],
                'desc': 'connect_calls_in_attach == total canonicalized direct attaches',
            },
            'total_connect_identity': {
                'lhs': total_connects, 'rhs': expected_total_connects,
                'pass': total_connects == expected_total_connects,
                'desc': 'total connects == canonicalized attaches + 2*bubbles',
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
        """Pretty-print the parallel-algorithm counters + identity check."""
        a = audit if audit is not None else self.audit

        lines = []
        title = "AUDIT 1 (parallel) -- discovery/commit counters"
        if header:
            title = f"{title}  ({header})"
        lines.append("=" * 64)
        lines.append(title)
        lines.append("=" * 64)
        lines.append("Discovery (requires instrument=True to be nonzero):")
        for k in ('discover_pruned', 'discover_attach_whole',
                  'discover_claim_empty', 'discover_demand_recorded'):
            lines.append(f"  {k:<28s} {a[k]:>14d}")
        lines.append("")
        lines.append("Commit (always on):")
        for k in ('nodes_with_demand_total', 'nodes_attached_total',
                  'nodes_covered_by_attach', 'extract_bubble_calls',
                  'mutations_moved', 'make_node_calls',
                  'connect_calls_in_attach', 'connect_calls_in_extract'):
            lines.append(f"  {k:<28s} {a[k]:>14d}")
        lines.append("")
        if a['recombine_calls']:
            lines.append("Per-recombine averages:")
            lines.append(f"  bubbles / offspring   {a['extract_bubble_calls']/a['recombine_calls']:>10.3f}")
            lines.append(f"  attaches / offspring  {a['nodes_attached_total']/a['recombine_calls']:>10.3f}")
        lines.append("")
        lines.append("Identities (via audit_check):")
        results = self.audit_check(audit=a, raise_on_fail=False)
        for r in results.values():
            mark = "OK" if r['pass'] else "FAIL"
            lines.append(f"  [{mark}] {r['desc']}: lhs={r['lhs']} rhs={r['rhs']}")
        lines.append("=" * 64)
        print("\n".join(lines))
