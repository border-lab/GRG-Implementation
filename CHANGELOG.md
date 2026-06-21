# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

For development workflow changes (testing, CI/CD, tooling), see [devtools/CHANGELOG.dev.md](devtools/CHANGELOG.dev.md).

## [Unreleased]

### Added

- `benchmark/benchmark_run.py` — single-run benchmark executor. Runs exactly one
  run of one method (GRG or NumPy) on one `.grg` file and writes a self-contained
  result JSON. Supports `--verification`, `--diagnostics`, `--serialize`, `--profile`,
  and `--seed` flags.
- `benchmark/benchmark_aggregate.py` — reads single-run JSONs produced by
  `benchmark_run.py`, groups by (file, method), computes summary statistics
  (mean/std/min/max), and writes combined CSV + JSON matching the original
  `benchmark_recombination.py` output format.
- `benchmark/haplotype_oracle.py` — per-offspring mutation correctness
  verification. For each offspring, walks up the GRG to collect inherited
  mutations and compares against the expected set derived from the parent
  haplotypes and recombination breakpoints. Integrated into `--verification`
  in both `benchmark_run.py` and `benchmark_recombination.py`. Also usable
  as a standalone CLI.

### Changed

- SLURM scripts (`sweep_tier1-4.sh`, `small_runs.sh`) now use per-run array jobs
  via `benchmark_run.py` instead of running all warmup + timed runs in a single
  monolithic process. Each array element is one (file, method, run_index) combination.
- Added companion `*_aggregate.sh` scripts for each tier, intended to run as
  `--dependency=afterok` jobs after the run arrays complete.
