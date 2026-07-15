# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

For development workflow changes (testing, CI/CD, tooling), see [devtools/CHANGELOG.dev.md](devtools/CHANGELOG.dev.md).

## [Unreleased]

### Added

- `benchmark/grg_recombination_parallel.py` — `ParallelNonDuplicationRecombination`,
  a node-aggregated discover+commit implementation of the parallel non-duplication
  recombination algorithm. A read-only, `ThreadPoolExecutor`-parallel discovery
  phase walks ancestral cones (pruned via per-node `Iu` spans) to record demanded
  genomic intervals and direct-attach edges per node, guarded by striped locks;
  a deterministic commit phase then splits old nodes and wires edges sequentially.
  Adds `benchmark/test_grg_recombination_parallel.py` (interval-algebra unit tests
  plus parallel-vs-sequential-vs-ground-truth equivalence tests against a mock
  `pygrgl`-shaped `FakeGRG`).
- `--method grg_parallel` / `--engine parallel` support across the benchmark
  tooling: `benchmark_run.py` (`--method grg_parallel`, `--max-workers`),
  `benchmark_recombination.py` (`--engine {sequential,parallel}`, `--max-workers`),
  and `multitree_check.py` (`--parallel`, `--max-workers`) all can now exercise
  `ParallelNonDuplicationRecombination` instead of the sequential reference.
  `benchmark_aggregate.py` recognizes the `grg_parallel` method and emits rows
  labeled `"GRG Parallel"` alongside `"GRG Native"`.
- `benchmark/benchmark_run.py` — single-run benchmark executor. Runs exactly one
  run of one method (GRG or NumPy) on one `.grg` file and writes a self-contained
  result JSON. Supports `--verification`, `--diagnostics`, `--serialize`, `--profile`,
  and `--seed` flags.
- `benchmark/benchmark_aggregate.py` — reads single-run JSONs produced by
  `benchmark_run.py`, groups by (file, method), computes summary statistics
  (mean/std/min/max), and writes combined CSV + JSON matching the original
  `benchmark_recombination.py` output format.

### Fixed

- `grg_recombination.py` / `grg_numpy_baseline.py` called `get_mutations_for_node()` /
  `get_node_mutation_pairs()` with an `allow_sort=False` keyword that is not part of
  upstream `pygrgl`'s API (verified against the official Python API docs and every
  local `grgl` source checkout available) and raised `TypeError` against a vanilla
  `pygrgl` install. Dropped the keyword at all four call sites; call-site behavior
  is unaffected since callers already re-sort the returned mutations by position in
  Python whenever order matters.
- Upstream `grgl` bug (`include/grgl/grg.h`, `MutableGRG::nodesAreTopo()`): unlike
  its sibling `nodesAreOrdered()`, `nodesAreTopo()` didn't account for `m_negativeNodes`
  being non-empty. Once any negative-ID node exists (i.e. after any
  `make_node(negative=True)`, which recombination calls once per offspring), the
  `save_grg`/`simplifyAndSerialize` fast path (`fastCompleteDFS` -> `getOrderedNodes()`)
  was taken instead of the safe generic DFS. `getOrderedNodes()` reconstructs order by
  assuming negative-node IDs are numerically contiguous at the end of ID space, which
  doesn't hold once negative and regular node IDs are interleaved across many
  recombination calls (the normal case: one negative offspring node then several
  regular bubble nodes, repeated per offspring). This silently duplicated some NodeIDs
  and dropped others from the reconstructed order, which then made the serializer visit
  a node twice and hit `release_assert(m_nodeIdMap[nodeId] == INVALID_NODE_ID)` — a hard
  process abort, reproducible with either recombination engine (sequential or parallel)
  after 2+ generations. Fixed by making `nodesAreTopo()` also require
  `m_negativeNodes.empty()`, matching `nodesAreOrdered()`; this trades the fast path for
  the always-correct generic DFS whenever negative nodes are present. Patched in the
  local `grgl` checkout and rebuilt/reinstalled into `.venv`'s `pygrgl`; not a change to
  this repo's own source (`grgl/` is a gitignored external dependency), but recorded
  here since it blocked multi-generation `--serialize`/post-run graph-size measurements
  in `benchmark_run.py` and `benchmark_recombination.py --serialize`.

### Changed

- SLURM scripts (`sweep_tier1-4.sh`, `small_runs.sh`) now use per-run array jobs
  via `benchmark_run.py` instead of running all warmup + timed runs in a single
  monolithic process. Each array element is one (file, method, run_index) combination.
- Added companion `*_aggregate.sh` scripts for each tier, intended to run as
  `--dependency=afterok` jobs after the run arrays complete.
