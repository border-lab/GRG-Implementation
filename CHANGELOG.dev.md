# Development Workflow Changelog

Changes to development tooling, CI/CD, testing infrastructure, and documentation systems.

## [Unreleased]

### Changed

- `simulate_grg_recombination` accepts optional `generation_index` and
  `offspring_ledger` parameters. When `offspring_ledger` is not None,
  per-offspring provenance (parents, segments, generation) is recorded
  for downstream verification. Zero overhead when unused.
- Benchmark pipeline restructured from monolithic to two-stage (run + aggregate).
  `benchmark_run.py` executes a single (file, method, run_index) job;
  `benchmark_aggregate.py` collects the per-run JSONs into the final CSV + JSON.
  This enables per-run SLURM parallelism, independent memory allocation per method,
  and fault tolerance (a failed run doesn't lose completed runs).
- SLURM sweep scripts updated to use array-job indexing over the
  (file, method, run) matrix with companion aggregation dependency jobs.
