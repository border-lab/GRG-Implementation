"""
grg_recombination_native.py

C++-backed drop-in replacement for grg_recombination.NonDuplicationRecombination.
Re-exports `simulate_grg_recombination` and `recombination_intervals` from the
Python module so existing callers see the same module-level API.

Usage:
    # In place of `from grg_recombination import NonDuplicationRecombination`:
    from grg_recombination_native import NonDuplicationRecombination
"""

from grg_recombination import (  # noqa: F401 -- re-exported for callers
    recombination_intervals,
    simulate_grg_recombination,
    compute_grg_structural_stats,
)

from grg_recomb_native._grg_recomb_native import (
    NonDuplicationRecombiner as _NativeRecombiner,
)


# Canonical key lists, used by `_fresh_audit` / `_fresh_stats` and as the
# source of truth for parity checks against the C++-backed audit/stats dicts.
# Order matches the Python NonDuplicationRecombination._fresh_audit() /
# _fresh_stats() definitions in grg_recombination.py.
_AUDIT_KEYS = (
    "pruning",
    "pruning_root",
    "path_compression",
    "path_compression_attach",
    "decomposition",
    "direct_attach",
    "direct_attach_root",
    "bubble_strip",
    "bubble_split",
    "direct_attach_dup",
    "bubble_fill",
    "bubble_strip_partial",
    "bubble_split_partial",
    "bubble_strip_partial_rt",
    "skip_empty_interval",
    "skip_already_visited",
    "skip_empty_trim",
    "visits",
    "extract_bubble_calls",
    "connect_calls_in_attach",
    "connect_calls_in_extract",
    "make_node_calls",
    "recombine_calls",
    "recurse_attach_calls",
)

_STATS_KEYS = (
    "init_caches_time",
    "recurse_attach_time",
    "apply_bubbles_time",
    "sync_to_grg_time",
    "clear_caches_time",
    "flush_samples_time",
    "sort_mutations_time",
    "get_mutation_by_id_time",
    "add_mutation_time",
    "remove_mutation_time",
    "make_node_time",
    "connect_time",
    "get_mutation_by_id_calls",
    "add_mutation_calls",
    "remove_mutation_calls",
    "make_node_calls",
    "connect_calls",
    "offspring_count",
    "recurse_attach_calls",
    "segments_processed",
    "visits_total",
    "bubbles_created",
    "mutations_moved",
)


class NonDuplicationRecombination:
    """Python wrapper over the C++ NonDuplicationRecombiner.

    Exposes the legacy Python NonDuplicationRecombination interface so that
    `simulate_grg_recombination` and `benchmark_recombination.py` are
    backend-agnostic. Maintains NEGATIVE_NODE_IDS Python-side (mirrored from
    the C++ side after each recombine) to avoid the O(N) per-access copy that
    routing through the C++ vector property would incur in the simulator's
    inner loop.
    """

    debug_mode = False  # class attribute for legacy compatibility

    def __init__(self, grg, instrument: bool = False):
        self._native = _NativeRecombiner(grg, instrument)
        # Mirror of self._native.negative_node_ids, grown incrementally
        # via self._native.last_raw_offspring_id after each recombine_multi.
        self.NEGATIVE_NODE_IDS = []

    # ------------------------------------------------------------------
    # Static factories matching Python NonDuplicationRecombination
    # ------------------------------------------------------------------

    @staticmethod
    def _fresh_audit():
        return {k: 0 for k in _AUDIT_KEYS}

    @staticmethod
    def _fresh_stats():
        # Mirror Python: timing fields are float (default 0.0), counters are int (0).
        return {k: (0.0 if k.endswith("_time") else 0) for k in _STATS_KEYS}

    # ------------------------------------------------------------------
    # Forwarded properties
    # ------------------------------------------------------------------

    @property
    def grg(self):
        return self._native.grg

    @property
    def audit(self):
        return self._native.audit  # Fresh dict each access (C++ -> py::dict).

    @property
    def stats(self):
        return self._native.stats

    @property
    def instrument(self):
        return self._native.instrument

    @property
    def defer_sample_updates(self):
        return self._native.defer_sample_updates

    @defer_sample_updates.setter
    def defer_sample_updates(self, value):
        self._native.defer_sample_updates = value

    # debug_mode is a class attribute (legacy); forward only when set on the
    # instance.
    def __setattr__(self, name, value):
        if name == "debug_mode":
            # Mirror onto the C++ side so the binding sees the same flag.
            object.__setattr__(self, name, value)
            if hasattr(self, "_native"):
                self._native.debug_mode = value
        else:
            object.__setattr__(self, name, value)

    # ------------------------------------------------------------------
    # Forwarded methods
    # ------------------------------------------------------------------

    def recombine(self, hap_a, hap_b, breakpoint_):
        offspring_id = self._native.recombine(hap_a, hap_b, breakpoint_)
        # Mirror the new offspring's raw NodeID onto the Python list.
        self.NEGATIVE_NODE_IDS.append(self._native.last_raw_offspring_id)
        return offspring_id

    def recombine_multi(self, segments):
        offspring_id = self._native.recombine_multi(segments)
        self.NEGATIVE_NODE_IDS.append(self._native.last_raw_offspring_id)
        return offspring_id

    def end_generation(self):
        self._native.end_generation()

    def clear_pending_sample_removals(self):
        self._native.clear_pending_sample_removals()

    def apply_pending_bubbles(self):
        self._native.apply_pending_bubbles()

    def flush_sample_updates(self):
        self._native.flush_sample_updates()

    def reset_audit(self):
        self._native.reset_audit()

    def audit_reset(self):
        # Python class uses this name; forward to the same C++ method.
        self._native.reset_audit()

    def reset_stats(self):
        self._native.reset_stats()

    # ------------------------------------------------------------------
    # Audit invariants / pretty-print (stay Python-side, operate on dicts)
    # ------------------------------------------------------------------

    def audit_check(self, audit=None, raise_on_fail=True):
        """Verify the three audit identities; same logic as the Python class."""
        a = audit if audit is not None else self.audit

        bubble_cases = (
            a["bubble_strip"]
            + a["bubble_split"]
            + a["bubble_fill"]
            + a["bubble_strip_partial"]
            + a["bubble_split_partial"]
            + a["bubble_strip_partial_rt"]
        )
        direct_firing = (
            a["direct_attach"] + a["direct_attach_root"] + a["path_compression_attach"]
        )
        total_connects = a["connect_calls_in_attach"] + a["connect_calls_in_extract"]
        expected_connects = 2 * a["extract_bubble_calls"] + direct_firing
        expected_make_nodes = a["recombine_calls"] + a["extract_bubble_calls"]

        results = {
            "bubble_identity": {
                "lhs": a["extract_bubble_calls"],
                "rhs": bubble_cases,
                "pass": a["extract_bubble_calls"] == bubble_cases,
                "desc": "extract_bubble_calls == sum of 6 bubble-case counters",
            },
            "connect_identity": {
                "lhs": total_connects,
                "rhs": expected_connects,
                "pass": total_connects == expected_connects,
                "desc": "total connect calls == 2*bubbles + firing direct attaches",
            },
            "make_node_identity": {
                "lhs": a["make_node_calls"],
                "rhs": expected_make_nodes,
                "pass": a["make_node_calls"] == expected_make_nodes,
                "desc": "make_node_calls == recombine_calls + bubbles",
            },
        }

        if raise_on_fail:
            for name, r in results.items():
                if not r["pass"]:
                    raise AssertionError(
                        f"audit_check FAIL [{name}]: {r['desc']} -- "
                        f"lhs={r['lhs']}, rhs={r['rhs']}, "
                        f"delta={r['lhs'] - r['rhs']}"
                    )
        return results

    def audit_summary(self, audit=None, header=None):
        """Pretty-print the case histogram + identity check.
        Delegates to the Python implementation by constructing a temporary
        Python NonDuplicationRecombination just for the formatter, since the
        formatter is purely string formatting and not perf-sensitive."""
        # Import lazily to avoid the cost when audit_summary isn't used.
        from grg_recombination import NonDuplicationRecombination as _PyImpl

        # _PyImpl.audit_summary doesn't actually need a constructed instance
        # for its core logic, but it is a method. Reuse its formatting code
        # by binding it to ourselves: it only reads `audit` arg + calls
        # audit_check, both of which work on us.
        return _PyImpl.audit_summary(self, audit=audit, header=header)
