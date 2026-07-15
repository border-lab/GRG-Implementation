# Development Workflow Changelog

Changes to development tooling, CI/CD, testing infrastructure, and documentation systems.

## [Unreleased]

### Changed

- Recreated repo-root `.venv` on **Python 3.9.25** (Homebrew `python@3.9`) to match
  the shared research freeze (was temporarily on 3.14). Install from
  `requirements-py39.txt` plus the local patched wheel in `grgl_wheels/`
  (`pygrgl-2.6-cp39-cp39-macosx_15_0_arm64.whl`, built from the sibling
  `../grg/grgl` checkout with the `nodesAreTopo` negative-node fix). The
  Akshay-path `pygrgl @ file:///Users/akshayanand/...` pin and the
  `grg_recomb_native` editable (`benchmark/native/`) were omitted — the latter
  directory is not present in this checkout.


- Benchmark pipeline restructured from monolithic to two-stage (run + aggregate).
  `benchmark_run.py` executes a single (file, method, run_index) job;
  `benchmark_aggregate.py` collects the per-run JSONs into the final CSV + JSON.
  This enables per-run SLURM parallelism, independent memory allocation per method,
  and fault tolerance (a failed run doesn't lose completed runs).
- SLURM sweep scripts updated to use array-job indexing over the
  (file, method, run) matrix with companion aggregation dependency jobs.
