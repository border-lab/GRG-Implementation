# Native C++ Port of `NonDuplicationRecombination` — Design

This document captures the design decisions for porting `benchmark/grg_recombination.py`'s
`NonDuplicationRecombination` class to a native C++ pybind11 extension. The motivation, scope,
and quantitative justification are in [`benchmark/recombination_time_breakdown.md`](../recombination_time_breakdown.md);
this document is the implementation contract.

The intent is to develop locally in this repo, then upstream into `aprilweilab/grgl` once the
API stabilizes and the perf numbers are settled. Files are organized and styled to match grgl's
conventions so the upstream migration is mostly `git mv`.

## Goals & non-goals

### In scope
- Port `_recurse_attach`, `_extract_bubble`, `_build_ancestral_caches`, `_apply_pending_bubbles`,
  `_clear_modified_caches`, and the per-node cache state to a single C++ class.
- Expose a pybind11 binding for the class with `recombine` / `recombine_multi` entry points.
- Preserve the existing 13-case audit counters so `audit_check` invariants still hold.
- Preserve the multitree-cardinality property (multitree_check.py must still pass).

### Out of scope (stays in Python)
- `recombination_intervals`, `simulate_grg_recombination`, breakpoint sampling — these are
  cheap glue (~3% of total runtime) and easier to iterate on in Python.
- `compute_grg_structural_stats`, `_summary` — diagnostic helpers, no perf pressure.
- `audit_summary` printer — purely formatting, operates on the audit dict.
- `benchmark_recombination.py` driver, `multitree_check.py` — unchanged.

## Hard constraints (inherited from grgl)

These determine what the C++ code may use and how it must be built. Violating any of these breaks
the upstream merge or the local build:

| Constraint | Source |
|---|---|
| `CMAKE_CXX_STANDARD 11` — no `std::optional`, `std::make_unique`, structured bindings, `if constexpr` | `grgl/CMakeLists.txt:77` |
| `-Werror` on all C++ | `grgl/CMakeLists.txt:74` |
| `clang-tidy` with `WarningsAsErrors: '*'` (large allowlist of disabled checks) | `grgl/.clang-tidy` |
| `clang-format --dry-run -Werror` must pass on all `.cpp`/`.h` | `grgl/format-check.sh` |
| Use bundled pybind11 at `grgl/third-party/pybind11/` (2.14.0.dev1) — pip-installed pybind11 will not share types | smoke test, this session |
| Compile defines: `-DCOMPACT_NODE_IDS`, `-DVCF_GZ_SUPPORT` | `grgl/CMakeLists.txt:72` |
| Linkage: against `libgrgl.a` + `libvbyte` (transitively pulled in by `csr_storage.h`) | this session |

### Style conventions

Pulled from `grgl/.clang-format` and `grgl/src/transform.cpp`:

- `ColumnLimit: 120`, `IndentWidth: 4`, no tabs
- `PointerAlignment: Left` (`MutableGRG* g`, `const MutableGRG& g`)
- `BreakBeforeBraces: Attach` (K&R)
- `NamespaceIndentation: None` — do not indent inside `namespace grgl { ... }`
- `BinPackArguments: false`, `BinPackParameters: false`
- `SortIncludes: CaseSensitive` — project headers first, std last
- Classes: `PascalCase`. Methods/free functions: `camelCase`. Members: `m_camelCase`.
  Shared-ptr aliases: `XxxPtr`. Header guards: `GRGL_FILENAME_H`. Files: `snake_case.cpp`/`.h`.
- Local constants: `SCREAMING_SNAKE` + `constexpr`
- Doxygen `@param`/`@return` on public declarations in headers; algorithmic prose comments above
  complex functions in `.cpp` files
- Prefer `static_cast` over C-style casts in new code (existing code has some C-style casts;
  do not introduce more)
- Use `release_assert(...)` for runtime invariants (matches existing code)

Run `clang-format -i` over all touched files before each commit. CI will reject the diff otherwise.

## Final file layout

Local development layout (this repo):

```
benchmark/native/
    DESIGN.md                       # this document
    CMakeLists.txt
    setup.py                        # invokes CMake, mirrors grgl/setup.py pattern
    pyproject.toml
    include/
        grgl/
            recombination.h         # public header — placed under grgl/ so the
                                    # upstream move is `git mv` with no #include edits
    src/
        recombination.cpp
        bindings.cpp                # PYBIND11_MODULE block
    test/
        test_recombination.cpp      # GoogleTest cases (local-only, not packaged)
```

The header lives at `include/grgl/recombination.h` rather than `include/recombination.h` so that
`#include "grgl/recombination.h"` works both locally and after upstreaming — no source edits
required when the files move.

Upstream layout after migration (`aprilweilab/grgl`):

```
grgl/include/grgl/recombination.h   # git mv from benchmark/native/include/grgl/
grgl/src/recombination.cpp          # git mv from benchmark/native/src/
grgl/src/python/_grgl.cpp           # extend with the binding block from bindings.cpp
grgl/test/unit/test_recombination.cpp  # git mv from benchmark/native/test/
grgl/CMakeLists.txt                 # add to GRGL_CORE_SOURCES + GRGL_PUBLIC_HEADERS + GRGL_TEST_SOURCES
```

## Class API

Single class `grgl::NonDuplicationRecombiner` in `namespace grgl`. Methods match the Python class
one-for-one, renamed to grgl's `camelCase` convention.

```cpp
// include/grgl/recombination.h
#ifndef GRGL_RECOMBINATION_H
#define GRGL_RECOMBINATION_H

#include "grg.h"

#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace grgl {

// NodeID, NodeIDSizeT come from grgl/common.h (uint32_t with -DCOMPACT_NODE_IDS).
// BpPosition, MutationId come from grgl/mutation.h (uint64_t and uint32_t respectively).
// We use grgl's existing typedefs throughout; no shadow aliases.

/**
 * Per-decision-case histogram. Cardinalities defined identically to the Python
 * audit dict in NonDuplicationRecombination._fresh_audit(); see audit_check()
 * for the three invariants these must satisfy.
 */
struct AuditCounters {
    uint64_t pruning = 0;
    uint64_t pruningRoot = 0;
    uint64_t pathCompression = 0;
    uint64_t pathCompressionAttach = 0;
    uint64_t decomposition = 0;
    uint64_t directAttach = 0;
    uint64_t directAttachRoot = 0;
    uint64_t bubbleStrip = 0;
    uint64_t bubbleSplit = 0;
    uint64_t directAttachDup = 0;
    uint64_t bubbleFill = 0;
    uint64_t bubbleStripPartial = 0;
    uint64_t bubbleSplitPartial = 0;
    uint64_t bubbleStripPartialRt = 0;
    uint64_t skipEmptyInterval = 0;
    uint64_t skipAlreadyVisited = 0;
    uint64_t skipEmptyTrim = 0;
    uint64_t visits = 0;
    uint64_t extractBubbleCalls = 0;
    uint64_t connectCallsInAttach = 0;
    uint64_t connectCallsInExtract = 0;
    uint64_t makeNodeCalls = 0;
    uint64_t recombineCalls = 0;
    uint64_t recurseAttachCalls = 0;
};

/**
 * Phase-level wallclock + per-C++-call timing, populated when the recombiner
 * was constructed with `instrument = true`. Field names match the Python
 * `stats` dict in NonDuplicationRecombination._fresh_stats() (camelCase here,
 * snake_case in the dict serialization for parity with the Python interface).
 *
 * Overhead when instrumentation is on: ~20-40 ns per std::chrono::steady_clock
 * probe (vs ~150 ns/probe for Python's time.perf_counter). With ~10 probes per
 * moved mutation, this is well under 1% of total runtime even at biobank scale.
 */
struct RecombinerStats {
    // Phase wallclock (seconds)
    double initCachesTime = 0.0;
    double recurseAttachTime = 0.0;
    double applyBubblesTime = 0.0;
    double syncToGrgTime = 0.0;
    double clearCachesTime = 0.0;
    double flushSamplesTime = 0.0;
    double sortMutationsTime = 0.0;
    // C++-call wallclock (seconds)
    double getMutationByIdTime = 0.0;
    double addMutationTime = 0.0;
    double removeMutationTime = 0.0;
    double makeNodeTime = 0.0;
    double connectTime = 0.0;
    // C++-call counts
    uint64_t getMutationByIdCalls = 0;
    uint64_t addMutationCalls = 0;
    uint64_t removeMutationCalls = 0;
    uint64_t makeNodeCalls = 0;
    uint64_t connectCalls = 0;
    // Algorithmic counters
    uint64_t offspringCount = 0;
    uint64_t recurseAttachCalls = 0;
    uint64_t segmentsProcessed = 0;
    uint64_t visitsTotal = 0;
    uint64_t bubblesCreated = 0;
    uint64_t mutationsMoved = 0;
};

/**
 * Non-duplication GRG recombination algorithm. Owns the MutableGRG it operates
 * on plus all per-node caches required by the iterative DFS. Constructing the
 * recombiner runs a one-pass O(V + E) initialization of all caches; subsequent
 * calls to recombine() / recombineMulti() execute entirely in C++ with no
 * Python callbacks (drop the GIL via py::gil_scoped_release in the binding).
 *
 * Invalidation rules (see recombination.cpp for full discussion):
 *   - extractBubble() invalidates node_id (lost mutations) and bubble_id (new slot)
 *   - extractBubble() does NOT invalidate node_id's parents (their span / anc_cov
 *     depend only on themselves and their own ancestors, both unchanged)
 *   - After sortMutations() at end-of-generation, mutationCache and posCache
 *     are cleared; position-derived caches (span, ancCov, upEdges) survive
 *     because positions do not change.
 */
class NonDuplicationRecombiner {
public:
    explicit NonDuplicationRecombiner(MutableGRGPtr grg, bool instrument = false);

    /** Two-parent recombine with a single crossover at `breakpoint`. */
    NodeID recombine(NodeID hapA, NodeID hapB, Position breakpoint);

    /**
     * Multi-segment recombine. `segments` is a list of (sourceParent, endCoord)
     * pairs covering [0, genomeLength); start coord is implicit from the
     * preceding segment's end.
     */
    NodeID recombineMulti(const std::vector<std::pair<NodeID, Position>>& segments);

    /** Apply pending bubble mutation moves; called automatically at end of recombine*. */
    void applyPendingBubbles();
    void flushSampleUpdates();

    /**
     * End-of-generation hook. Encapsulates the post-loop bookkeeping that
     * `simulate_grg_recombination` currently does in Python:
     *   1. grg->sortMutations() (compacts soft-deleted mutation entries)
     *   2. Clear mutationCache + posCache (sortMutations renumbers MutationIds)
     *   3. Record sortMutationsTime if instrumented
     * The simulator must call this once per generation. Position-derived caches
     * (span, ancCov, upEdges) intentionally survive because positions don't
     * change across sortMutations.
     */
    void endGeneration();

    /** Clear pending sample removals (called after wholesale set_samples). */
    void clearPendingSampleRemovals();

    void setDeferSampleUpdates(bool defer) { m_deferSampleUpdates = defer; }
    bool getDeferSampleUpdates() const { return m_deferSampleUpdates; }

    void setDebugMode(bool d) { m_debug = d; }
    bool getDebugMode() const { return m_debug; }

    bool isInstrumented() const { return m_instrument; }

    const AuditCounters& getAudit() const { return m_audit; }
    void resetAudit() { m_audit = AuditCounters{}; }

    const RecombinerStats& getStats() const { return m_stats; }
    /** Zero both stats and audit (matches Python `reset_stats()`). */
    void resetStats() { m_stats = RecombinerStats{}; m_audit = AuditCounters{}; }

    /** Match Python NEGATIVE_NODE_IDS — see "Offspring ID convention" below. */
    const std::vector<NodeID>& getNegativeNodeIds() const { return m_negativeNodeIds; }

private:
    // Construction
    void buildAncestralCaches();
    void growNodeArrays(NodeID nodeId);
    void syncToGrg();

    // Hot path
    void recurseAttach(NodeID root, NodeID offspring, Position L, Position R);
    NodeID extractBubble(NodeID nodeId,
                         const std::vector<MutationId>& relMutIds,
                         NodeID offspringId);

    // Cache access
    const std::vector<std::pair<MutationId, Position>>& getNodeMutations(NodeID nodeId);
    const std::vector<NodeID>& getUpEdgesCached(NodeID nodeId);
    IntervalOpt getNodeAndAncestorSpan(NodeID nodeId);
    IntervalOpt getAncestralCoverage(NodeID nodeId);

    // Cleanup
    void clearModifiedCaches();
    NodeID registerOffspring(NodeID offspringId);

    MutableGRGPtr m_grg;
    Position m_genomeLength;

    // Per-node arrays (dense, indexed by NodeID, grown as bubbles/offspring added)
    std::vector<IntervalOpt> m_spanCache;
    std::vector<IntervalOpt> m_ancCovCache;
    std::vector<uint64_t> m_visitedGen;
    std::vector<uint64_t> m_connectedGen;
    uint64_t m_genVisited = 0;
    uint64_t m_genConnected = 0;

    // Hash-keyed caches (eligible nodes may not yet have IDs at init time)
    std::unordered_map<NodeID, std::vector<std::pair<MutationId, Position>>> m_mutationCache;
    std::unordered_map<NodeID, std::vector<Position>> m_posCache;
    std::unordered_map<NodeID, std::vector<NodeID>> m_upEdgesCache;

    // Deferred work
    struct BubbleOp {
        NodeID nodeId;
        NodeID bubbleId;
        std::vector<MutationId> relMutIds;
    };
    std::vector<BubbleOp> m_pendingBubbles;
    std::unordered_set<NodeID> m_pendingSampleRemovals;
    std::unordered_set<NodeID> m_modifiedNodes;

    AuditCounters m_audit;
    RecombinerStats m_stats;
    bool m_instrument = false;
    bool m_debug = false;
    bool m_deferSampleUpdates = false;

    // Offspring tracking (mirror of Python NEGATIVE_NODE_IDS + _negative_node_index)
    std::vector<NodeID> m_negativeNodeIds;
    std::unordered_map<NodeID, size_t> m_negativeNodeIndex;
};

}  // namespace grgl

#endif  // GRGL_RECOMBINATION_H
```

`IntervalOpt` is defined in the implementation; see "Cache representation" below.

## Cache representation — `IntervalOpt`

Python's three-state for `span_cache[node]` and `anc_cov_cache[node]` (`False` = uninitialized,
`None` = no parents, `(lo, hi)` = computed) needs a C++11 equivalent. Decision: small POD struct,
declared in the header.

```cpp
struct IntervalOpt {
    Position lo;
    Position hi;
    enum class State : uint8_t { Uninit = 0, Empty = 1, Set = 2 };
    State state;

    static IntervalOpt uninit() { return IntervalOpt{0, 0, State::Uninit}; }
    static IntervalOpt empty()  { return IntervalOpt{0, 0, State::Empty}; }
    static IntervalOpt set(Position lo, Position hi) { return IntervalOpt{lo, hi, State::Set}; }

    bool isUninit() const { return state == State::Uninit; }
    bool isEmpty()  const { return state == State::Empty; }
    bool isSet()    const { return state == State::Set; }
};
```

Trade-off vs. parallel-vector flag arrays: the struct is 24 bytes per node (with padding) vs.
17 bytes for two `uint64_t` + one `uint8_t` flag spread across three arrays. At ~2.5M nodes this
is ~17 MB extra — negligible against the ~700 MB working set the largest benchmark already
uses. Take the readability.

## NodeID-keyed cache structure

Decision: mirror the Python layout exactly for the first cut.

| Cache | Type | Why |
|---|---|---|
| `m_spanCache`, `m_ancCovCache` | `std::vector<IntervalOpt>` | All input-graph nodes are densely populated by the init pass; bubbles/offspring extend the vector via `growNodeArrays` |
| `m_visitedGen`, `m_connectedGen` | `std::vector<uint64_t>` | Same |
| `m_mutationCache`, `m_posCache`, `m_upEdgesCache` | `std::unordered_map<NodeID, ...>` | Eligible for eviction by `clearModifiedCaches`; sparse re-population after `sortMutations` |
| `m_pendingSampleRemovals`, `m_modifiedNodes` | `std::unordered_set<NodeID>` | Small per-recombine, mirror Python sets |

Flattening the hash-keyed caches into struct-of-arrays for cache locality is a known follow-up
(deferred optimization #2 in the time-breakdown report). Land correctness first, profile, then
revisit.

## Offspring ID convention

Decision: **C++ returns the raw `NodeID` (negative); Python wrapper maintains the
`NEGATIVE_NODE_IDS` reverse-lookup table.**

Rationale: the negative-ID table is purely a Python-side bookkeeping concern (the simulator uses
it for stable offspring identifiers across generations). Keeping it in Python avoids exposing
a stateful table through the binding and makes the C++ class re-entrant.

The C++ class still maintains `m_negativeNodeIds` because `_register_offspring` is called from
inside the recombine path in the current Python code; for parity the C++ side mirrors it and
exposes via `getNegativeNodeIds()`. The Python wrapper consumes that for `audit_summary` and
ID translation.

## Audit / instrumentation exposure

Decision: **expose `AuditCounters` to Python as a `dict` via a `to_dict()` method on the binding.**

`audit_summary` and `audit_check` in Python operate on a `dict[str, int]`. By returning a dict
from the binding with field names matching the existing Python audit keys (`pruning`,
`pruning_root`, etc. — snake_case in the dict even though C++ field names are camelCase), we keep
the printing/identity-check helpers byte-for-byte unchanged.

Field-name mapping in `bindings.cpp`:

```cpp
py::dict toDict(const AuditCounters& a) {
    py::dict d;
    d["pruning"]                    = a.pruning;
    d["pruning_root"]               = a.pruningRoot;
    d["path_compression"]           = a.pathCompression;
    d["path_compression_attach"]    = a.pathCompressionAttach;
    d["decomposition"]              = a.decomposition;
    d["direct_attach"]              = a.directAttach;
    d["direct_attach_root"]         = a.directAttachRoot;
    d["bubble_strip"]               = a.bubbleStrip;
    d["bubble_split"]               = a.bubbleSplit;
    d["direct_attach_dup"]          = a.directAttachDup;
    d["bubble_fill"]                = a.bubbleFill;
    d["bubble_strip_partial"]       = a.bubbleStripPartial;
    d["bubble_split_partial"]       = a.bubbleSplitPartial;
    d["bubble_strip_partial_rt"]    = a.bubbleStripPartialRt;
    d["skip_empty_interval"]        = a.skipEmptyInterval;
    d["skip_already_visited"]       = a.skipAlreadyVisited;
    d["skip_empty_trim"]            = a.skipEmptyTrim;
    d["visits"]                     = a.visits;
    d["extract_bubble_calls"]       = a.extractBubbleCalls;
    d["connect_calls_in_attach"]    = a.connectCallsInAttach;
    d["connect_calls_in_extract"]   = a.connectCallsInExtract;
    d["make_node_calls"]            = a.makeNodeCalls;
    d["recombine_calls"]            = a.recombineCalls;
    d["recurse_attach_calls"]       = a.recurseAttachCalls;
    return d;
}
```

A `RecombinerStats` → `py::dict` helper follows the same pattern, with snake_case keys matching
the existing Python `stats` keys (`init_caches_time`, `recurse_attach_time`, `apply_bubbles_time`,
`sync_to_grg_time`, `clear_caches_time`, `flush_samples_time`, `sort_mutations_time`,
`get_mutation_by_id_time`, `add_mutation_time`, `remove_mutation_time`, `make_node_time`,
`connect_time`, `get_mutation_by_id_calls`, `add_mutation_calls`, `remove_mutation_calls`,
`make_node_calls`, `connect_calls`, `offspring_count`, `recurse_attach_calls`,
`segments_processed`, `visits_total`, `bubbles_created`, `mutations_moved`).

## Public API parity matrix

The benchmark stack (`benchmark_recombination.py` + `simulate_grg_recombination` +
`multitree_check.py`) must observe identical behavior whether routed through the Python
`NonDuplicationRecombination` or the C++ port wrapped in a Python shim. The full surface area
the benchmark stack touches, and where each lives in the new layout:

### Constructor

| Call | Wrapper does |
|---|---|
| `NonDuplicationRecombination(g)` | `__init__(self, g)` → constructs C++ `NonDuplicationRecombiner(g, False)` |
| `NonDuplicationRecombination(g, instrument=True)` | `__init__(self, g, instrument=True)` → constructs C++ `NonDuplicationRecombiner(g, True)` |

### Attributes read directly from the recombiner

| Python access pattern | Source (line) | Backed by |
|---|---|---|
| `recomb.grg` | `simulate_grg_recombination:1067, 1093, 1102` | Python wrapper holds the `pygrgl.MutableGRG` reference; forwards |
| `recomb.defer_sample_updates` (read/write) | `simulate_grg_recombination:1075-76, 1113` | Wrapper property → C++ `getDeferSampleUpdates` / `setDeferSampleUpdates` |
| `recomb.instrument` | `simulate_grg_recombination:1100, 1103` | Wrapper property → C++ `isInstrumented()` |
| `recomb.NEGATIVE_NODE_IDS` | `simulate_grg_recombination:1089` | Wrapper attribute (Python-side, per offspring-ID decision); mirrors C++ `getNegativeNodeIds()` |
| `recomb.audit` | `benchmark_recombination.py:609`, `audit_check` impl | Wrapper property → C++ `getAudit()` → `audit_dict()` |
| `recomb.stats` | `benchmark_recombination.py:594, 605`, `simulate_grg_recombination:1104` | Wrapper property → C++ `getStats()` → `stats_dict()`. Must support both `recomb.stats["k"]` access and `dict(recomb.stats)` snapshot |

### Methods called on the recombiner

| Python call | Backed by |
|---|---|
| `recomb.recombine_multi(segments)` | C++ `recombineMulti(segments)` (with GIL released) |
| `recomb.audit_check(raise_on_fail=False)` | Stays in Python wrapper; reads `recomb.audit` dict and checks the 3 invariants |
| `recomb.audit_summary(audit=..., header=...)` | Stays in Python wrapper; pretty-print operating on a dict |
| `recomb.reset_stats()` | C++ `resetStats()` (zeros both `m_stats` and `m_audit`) |

### Static methods on the class

| Python call | Backed by |
|---|---|
| `NonDuplicationRecombination._fresh_audit()` | Wrapper `@staticmethod` returning a zero-initialized dict with all audit keys (used for aggregation in `--diagnostics`) |
| `NonDuplicationRecombination._fresh_stats()` | Wrapper `@staticmethod` returning a zero-initialized dict with all stats keys |

### Private state touched directly by `simulate_grg_recombination` (current code)

These currently reach into class privates. They need either equivalent C++ methods or a wrapper
refactor that calls a single high-level method. **Decision: collapse all four into an
`endGeneration()` method on the C++ class plus an explicit `clearPendingSampleRemovals()`.**

| Current direct access | Replaced by |
|---|---|
| `recomb._pending_sample_removals.clear()` (line 1095) | C++ `clearPendingSampleRemovals()` |
| `recomb.grg.sort_mutations()` (line 1102) | C++ `endGeneration()` (internal call) |
| `recomb._mutation_cache.clear()` (line 1110) | C++ `endGeneration()` (internal call) |
| `recomb._pos_cache.clear()` (line 1111) | C++ `endGeneration()` (internal call) |
| `recomb.stats["sort_mutations_time"] += ...` (line 1104) | C++ `endGeneration()` records the timing internally when instrumented |

After the refactor, `simulate_grg_recombination` is changed to call
`recomb.clear_pending_sample_removals()` (after `set_samples`) and `recomb.end_generation()`
(in place of the manual `sort_mutations + cache clear + stats increment` block). Both methods
exist on both the Python and the C++ wrapper, so the simulator is implementation-agnostic.

### Methods on `recomb.grg` used by the simulator

Forwarded through the wrapper's `grg` attribute (already on the `pygrgl.MutableGRG` API):
`get_sample_nodes()`, `set_samples(...)`, `num_nodes`, `num_edges`. No new bindings needed.

### Benchmark-driver flags vs class-internal flags

Flags that touch class behavior (must port) vs orchestration flags (Python-only):

| Flag | Class-internal? | Routing |
|---|---|---|
| `--diagnostics` | YES — requires `instrument=True` on the recombiner | C++ instrumentation must be live; wrapper exposes `stats`/`audit` dicts |
| `--verification` | YES — calls `audit_check` | Wrapper's `audit_check` reads C++ audit dict |
| `--profile` | NO — wraps the simulator in cProfile | Unchanged |
| `--serialize` | NO — pygrgl save/load on the GRG | Unchanged |
| `--skip-numpy` | NO | Unchanged |
| `--warmup`, `--runs`, `--num-generations`, `--offspring-per-couple` | NO | Unchanged |

## Instrumentation (`instrument=True` path)

Decision: **port the full instrumentation surface to C++. `--diagnostics` must produce the same
JSON shape regardless of which implementation backs the recombiner.**

The benchmark driver `benchmark_recombination.py` reads from `recomb.stats[<key>]` at multiple
points (line refs in the current code):

- `diag_recomb.stats["init_caches_time"]` (line 594) — captured once at ctor
- `dict(diag_recomb.stats)` per generation (line 605) — full snapshot for JSON output
- `snapshot["remove_mutation_time"] / snapshot["remove_mutation_calls"]` (lines 629–633) — per-call cost prints
- `snapshot["add_mutation_time"] / snapshot["add_mutation_calls"]` (lines 634–636)
- `snapshot["offspring_count"]`, `["bubbles_created"]`, `["mutations_moved"]`, `["visits_total"]` (lines 624–628)

Dropping `stats` would force a parallel rewrite of `--diagnostics` and break the JSON contract
with downstream tooling.

Justification for keeping despite the time-breakdown report concluding the wallclock breakdowns
weren't needed for the port itself:

1. **`std::chrono::steady_clock::now()` is cheap** — ~20-40 ns on Apple Silicon vs ~150 ns for
   Python's `time.perf_counter`. With ~10 probes per moved mutation, instrumentation overhead is
   <1% of total runtime even when on. The Python-side cost (3-7% of runtime per the report) was
   what made dropping it tempting in Python; that cost mostly disappears in C++.
2. **`--diagnostics` is the load-bearing tool for future tuning.** Next bottleneck (likely
   `MutableGRG` mutator calls per the time-breakdown analysis) needs the per-call cost breakdown
   to diagnose.
3. **Cost to keep is bounded** — ~50 lines of `if (m_instrument) { ... }` blocks in
   `recombination.cpp`, matching the same probe points the Python code has. The plumbing already
   exists conceptually.

Exposure: `getStats()` returns the C++ struct; `bindings.cpp` provides a `stats_dict()` method
that returns a `py::dict` with snake_case keys matching the existing Python `stats` keys. Python
wrapper exposes both `recomb.stats` (calling `stats_dict()` on access) and `dict(recomb.stats)`
semantics via `__getitem__`. See "Public API parity matrix" below for the full key list.

## GIL handling

Decision: release the GIL around the C++ entry points in `bindings.cpp`.

```cpp
.def("recombine_multi",
     [](NonDuplicationRecombiner& self,
        const std::vector<std::pair<NodeID, Position>>& segments) {
         py::gil_scoped_release release;
         return self.recombineMulti(segments);
     })
```

No Python callbacks fire inside the inner loop (all data is C++-owned post-construction), so the
GIL can be released for the full duration. This is needed for the routine to compose with
multi-threaded Python drivers later; it does not by itself parallelize a single recombine call.

## Build setup

Decision: **hand-rolled `setup.py` that invokes CMake, mirroring `grgl/setup.py`.**

This is the same pattern grgl uses (a `CMakeExtension` subclass + `build_ext` override). Mirroring
it now means the upstream merge moves the binding registration into `_grgl.cpp` and deletes our
`setup.py` / `CMakeLists.txt` — nothing structurally new.

CMake responsibilities:
1. Pull in grgl's source tree via `add_subdirectory(${GRGL_ROOT})` where `GRGL_ROOT` defaults to
   `../../grgl` and is overridable via `-DGRGL_ROOT=...`. This builds `libgrgl.a` + `libvbyte`
   as side effects, ensuring our extension links against the same compiled grgl that pygrgl
   does (no ABI drift).
2. Use grgl's bundled pybind11 at `${GRGL_ROOT}/third-party/pybind11`.
3. Apply `-DCOMPACT_NODE_IDS`, `-Werror`, `-std=c++11` to match grgl's compile flags.
4. Build target `_grg_recomb_native` as a pybind11 module, link against `grgl` + `vbyte`.

`setup.py` responsibilities:
1. `pip install -e .` triggers a CMake build into a temp directory, then copies the `.so` into
   the package install location.
2. Mirror `grgl/setup.py`'s `CMakeExtension` + `CMakeBuild` pattern.

Python package name: `grg_recomb_native`. From Python:

```python
import pygrgl
from grg_recomb_native import NonDuplicationRecombiner

grg = pygrgl.load_mutable_grg("example.grg")
recomb = NonDuplicationRecombiner(grg)
offspring_id = recomb.recombine_multi(segments)
```

The Python wrapper in `benchmark/grg_recombination.py` becomes a thin compatibility shim that
holds a `grg_recomb_native.NonDuplicationRecombiner` and exposes `audit_summary` / `audit_check`
(operating on the dict returned by `recomb.audit`) plus the `NEGATIVE_NODE_IDS` reverse-lookup.
Existing call sites in `benchmark_recombination.py` should not need changes.

## Test plan

### C++ unit tests (`test/test_recombination.cpp`)

GoogleTest, locally-run only (not built when the package is installed via pip):

1. **Audit identity #1**: `extractBubbleCalls == sum of 6 bubble-case counters`.
2. **Audit identity #2**: `connectCallsInAttach + connectCallsInExtract == 2 * extractBubbleCalls + (directAttach + directAttachRoot + pathCompressionAttach)`.
3. **Audit identity #3**: `makeNodeCalls == recombineCalls + extractBubbleCalls`.
4. **Fixed-graph regression**: hand-built 8-node GRG (matching `test.py` at repo root), recombine
   with a known segment list, assert specific offspring topology.
5. **One-pass cache init invariant**: post-construction, `m_spanCache[node].state == Set` (or
   `Empty` for nodes with no mutations and no parent spans) for every input-graph node.

### Python parity tests (`benchmark/native/test_parity.py`)

For a small set of `.grg` files (e.g. `50inds_1k_snps.grg`, `500inds_10k_snps.grg`), build both
the Python `NonDuplicationRecombination` and the C++ `NonDuplicationRecombiner` from the same
graph; recombine with a fixed PRNG seed; assert:

1. Identical offspring NodeID assignments.
2. Identical audit dict (every key, every count).
3. Identical post-recombination GRG state (num_nodes, num_edges, mutation counts).
4. Identical `stats` dict keys (counts must match exactly when both run with `instrument=True`;
   timing fields will differ in absolute value but their existence and structure must match).
5. `_fresh_audit()` and `_fresh_stats()` return identical key sets on both implementations.

This is the strongest regression net we get — the audit identities only check that the C++ code
is internally consistent; parity tests confirm it matches the validated Python reference.

### `--diagnostics` JSON parity test

Run `benchmark_recombination.py --diagnostics --num-generations 1` against the same `.grg` with
both implementations; assert the resulting JSON files have:

- Identical `diagnostics.<file>.per_generation` key sets (every snapshot field).
- Identical `audit_aggregated` dicts (full key/value match).
- `init_caches_time_s` present on both (values may differ).

This protects the downstream JSON contract — anything consuming the benchmark output should not
need to know which implementation generated it.

### End-to-end

`benchmark_recombination.py` runs unchanged against the wrapper. The `--verification` flag
exercises audit identities + multitree-cardinality check; these must pass on the C++ path with
no audit-counter delta vs. the Python path. The `--diagnostics` flag exercises the full `stats`
plumbing; must produce a JSON-comparable output to the Python path.

## Verification before declaring done

1. `clang-format --dry-run -Werror` clean on all `.cpp`/`.h`.
2. `clang-tidy` clean (run via grgl's `format-check.sh` after upstream move; locally, can run
   manually).
3. Parity tests (offspring IDs, audit, stats keys, fresh-dict helpers) pass on the three smallest
   `.grg` files in `grg_files/`.
4. `--verification` passes on at least one biobank-scale file (`32000inds_1m_snps.grg`).
5. `--diagnostics` JSON parity holds against the Python implementation on at least one small
   `.grg` (sub-second runtime so iteration is cheap).
6. `mean_time_ms` improvement on `64000inds_1m_snps.grg` ≥ 3× vs. the Python reference, per the
   time-breakdown report's prediction band.
7. Multitree-cardinality check passes (0 violations).

Once all seven are green, the upstream PR is essentially preparing the move.

## Migration path to grgl upstream

Steps when ready to upstream:

1. `git mv benchmark/native/include/grgl/recombination.h grgl/include/grgl/recombination.h`
2. `git mv benchmark/native/src/recombination.cpp grgl/src/recombination.cpp`
3. Fold `benchmark/native/src/bindings.cpp`'s `py::class_<NonDuplicationRecombiner>` block into
   `grgl/src/python/_grgl.cpp` (after the existing `MutableGRG` binding).
4. `git mv benchmark/native/test/test_recombination.cpp grgl/test/unit/test_recombination.cpp`
5. Edit `grgl/CMakeLists.txt`:
   - Append `src/recombination.cpp` to `GRGL_CORE_SOURCES`
   - Append `include/grgl/recombination.h` to `GRGL_PUBLIC_HEADERS`
   - Append `test/unit/test_recombination.cpp` to `GRGL_TEST_SOURCES`
6. Delete `benchmark/native/`, replace pure-Python `NonDuplicationRecombination` with shim
   importing from `pygrgl` (now exporting `NonDuplicationRecombiner`).
7. Bump `grgl/include/grgl/version.h` to `2.7` (minor — additive change).

No logic changes during the move; no `#include` edits (the header path is identical pre- and
post-move because we placed it under `include/grgl/` from day one).

## Open questions deferred (not blocking start of work)

These are flagged in the time-breakdown report or implied by the design, but do not need
resolution before code starts:

- **Cache flattening to SoA**: hash-keyed mutation/pos/up-edges caches → flat vectors with
  per-node offsets. Profile-driven, after correctness.
- **Bulk recombine API**: a single `recombineBatch(parents, segmentsPerOffspring)` entry point
  that holds the GIL release across many offspring. Would amortize the per-call binding overhead
  (~140 ns measured) if the batch loop turns out to matter.
- **`uint64_t` generation counter overflow**: the gen counter increments per recombine call;
  even at 10^9 recombines per process lifetime it does not wrap. No mitigation planned.
- **`sortMutations()` re-cache strategy**: currently clears `_mutation_cache` + `_pos_cache`
  wholesale at end-of-generation. Could be smarter (only clear nodes whose mutations actually
  moved), but the wholesale clear has not shown up as a bottleneck.
