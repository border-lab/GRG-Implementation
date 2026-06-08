# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project context

Research implementation around **Genotype Representation Graphs (GRG)**, focused on simulating recombination directly on the graph without materializing dense matrices.

**Background PDFs at the repo root** (read these first for context):
- `GRG Paper.pdf` — the foundational GRG paper. Explains what a GRG is, how it represents genotype data, and the construction algorithms.
- `GRG_Recombination-2.pdf` — describes the current recombination algorithm implemented in this repo, specifically the **Insertion method** (the "Non-duplication" / bubble-extraction approach realised by `NonDuplicationRecombination` in `benchmark/grg_recombination.py`).

There are **two parallel implementations** that share concepts but not code:

- **Repo-root pure-Python reference** (`grg.py`, `recombination.py`, `dot_product.py`, `plot_grg.py`): a small, readable `Node` / `GRG` model using dict-based parent/child adjacency with `[L, R)` interval edges. Used for the toy datasets in `test.py` / `test_toy_data.py` / `toy_data.json` and exploratory notebooks. **Not** what the benchmarks measure.
- **`benchmark/` production-grade implementation** built on top of `pygrgl` (the upstream C++ library exposed as Python). `grg_recombination.NonDuplicationRecombination` is the real algorithm — iterative DFS in `_recurse_attach`, per-node mutation / span / ancestral-coverage caches, deferred bubble application, and a one-pass O(V+E) Kahn's-sort cache init. This is what the paper figures and `.md` reports describe.

The `grgl/` subdirectory is an embedded clone of the upstream [aprilweilab/grgl](https://github.com/aprilweilab/grgl) C++/Python library. It is gitignored — treat it as read-only reference. `pygrgl` is installed into `.venv` separately.

## Environment

**Always activate `.venv` (Python 3.9) before running anything that imports `pygrgl`, msprime, or anything under `benchmark/`.** The venv lives at the repo root.

```bash
source .venv/bin/activate
```

Or invoke the interpreter directly: `.venv/bin/python <script>` (from `benchmark/`, use `../.venv/bin/python`).

`benchmark/requirements.txt` lists *only* the minimal extras (`numpy`, `tqdm`, `packaging`); `pygrgl`, `msprime`, `tskit`, `psutil`, `joblib` are already installed in `.venv`.

## Common commands

### Toy / repo-root scripts
```bash
.venv/bin/python test.py           # 8-node reference GRG; runs recombine() and plots
.venv/bin/python test_toy_data.py  # builds a GRG from toy_data.json; matrix + X^T v
```

### Generate `.grg` test inputs (msprime → tskit → grg convert/construct)
```bash
cd benchmark
../.venv/bin/python generate_grg.py -n 4000 -m 500000 --method ts
# emits ./grg_files/4000inds_500k_snps.grg
```
Key flags: `--input-trees <file.trees>` to skip msprime; `--method {ts,vcf,trees}` picks the conversion path; segment-size auto-derived as `max(200_000, 12.5 * num_samples)` (see `_derive_segment_size`).

### Run the recombination benchmark
```bash
cd benchmark
../.venv/bin/python benchmark_recombination.py \
    --grg-dir ../grg_files \
    --warmup 1 --runs 3 \
    --num-generations 2 \
    --offspring-per-couple 2 \
    --output-dir ./output
```
Outputs CSV + JSON to `--output-dir` (default `.`). Off-the-timed-path diagnostic flags (each adds wallclock but **does not** affect the headline mean_time_ms):

- `--skip-numpy` — skip the dense NumPy baseline (memory blows up past ~16k inds at 1M SNPs).
- `--diagnostics` — instrumented diagnostic pass (phase breakdowns, per-C++-call costs, audit-1 histogram), structural-stats fingerprint, and per-generation liveness/deadweight snapshots.
- `--verification` — after each non-warmup run: (1) audit identity checks on accumulated counters; (2) per-offspring multitree cardinality check (Approach B in `multitree_check.py`). Expensive at biobank scale.
- `--serialize` — after each generation of the last timed run, do `pygrgl.save_grg + load_mutable_grg`, report (nodes, edges) pre/post and save+load wallclock. Passive (graph is **not** replaced with the reloaded version).
- `--profile` — one extra generation under `cProfile`, top 30 by cumtime + tottime.

### Smoke-test a single file
```bash
cd benchmark
../.venv/bin/python benchmark_recombination.py \
    --grg-dir ../grg_files/50inds_1k_snps.grg \
    --warmup 0 --runs 1 --num-generations 1 --skip-numpy
```

### Standalone multitree check
```bash
cd benchmark
../.venv/bin/python multitree_check.py ../grg_files/<file>.grg --offspring-per-couple 2
```

There are no unit tests — the test surface is the audit-identity checks + multitree cardinality check exercised via `--verification`, plus `manual_toy_data_test.py` and the `test*.py` scripts at the root for the pure-Python reference.

## Architecture notes for the benchmark recombination

These are the load-bearing invariants. Touching them silently breaks correctness.

**13-case decision matrix in `_recurse_attach`** (`benchmark/grg_recombination.py:545`). Every visit past the skip guards lands in exactly one of 13 audit buckets (9 standard + 4 root variants). The matrix axes:
- *Relevant-mutation status*: `has_no_relevant`, `has_partial_relevant`, `has_all_relevant`.
- *Ancestral coverage of `[L, R)` relative to the node's `Iu = (anc_min, anc_max+1)`*: `full_coverage`, `ancestral_disjoint`, partial overlap, or `Iu is None` (root).

**Three audit identities** (`audit_check`, line 931):
1. `extract_bubble_calls == sum of 6 bubble-case counters` — every bubble was classified.
2. `total_connects == 2*bubbles + firing_direct_attaches` — `connect()` only fires in the prescribed paths.
3. `make_node_calls == recombine_calls + extract_bubble_calls` — no stray nodes.

`direct_attach_dup` (skipped by the `connected_gen` guard) is excluded from #2; `path_compression_attach` joins the direct-firing term because it fires one `connect()`.

**Caches and what invalidates them**:
- `_up_edges_cache`, `_mutation_cache`, `_pos_cache`, `span_cache`, `anc_cov_cache` are all populated up-front by `_build_ancestral_caches` via Kahn's sort.
- `_extract_bubble` only invalidates `node_id` (lost mutations) and `bubble_id` (fresh slot), **not** `node_id`'s parents (their span / anc_cov depend only on themselves and their own ancestors).
- After `grg.sort_mutations()` at end-of-generation, only `_mutation_cache` + `_pos_cache` are cleared; position-derived caches survive because positions don't change.

**Sample-set bookkeeping**: `simulate_grg_recombination` sets `defer_sample_updates=True`, then calls `grg.set_samples(new_offspring_ids)` once at the end. Per-offspring `_pending_sample_removals` is then dropped wholesale. Do not call `set_samples` inside the per-offspring loop.

**Negative node IDs**: offspring use negative IDs by convention (`grg.make_node(negative=True)`). `NEGATIVE_NODE_IDS` + `_negative_node_index` give the O(1) reverse map; `_register_offspring` returns the negative external handle.

## Conventions

- `.grg` files in `grg_files/` follow `<N>inds_<M>_snps.grg` so `_parse_expected_size` in the benchmarker can extract expected sizes. `k`/`m` suffixes (`500k`, `1m`) are honored.
- Benchmark `.md` reports (`benchmark/*_comparison.md`, `*_report.md`, etc.) are findings from prior runs. They are gitignored — treat them as historical context, not source of truth.
- `recombination*.ipynb` notebooks are exploratory and predate the current `grg_recombination.py`; assume they are out of date unless the user says otherwise.
- `pseudocode.py` is human-readable pseudocode (not runnable) for the original recurse-attach design.

## Things that bite

- Don't run benchmark scripts without sourcing `.venv` first — `pygrgl` will be missing and msprime/tskit will resolve to the wrong Python.
- The repo-root `grg.py` is **not** what `benchmark/grg_recombination.py` uses; don't conflate the two when editing.
- `grg_files/` and `benchmark/output/` are gitignored; expect them to be empty on a fresh clone (use `generate_grg.py`).
- The diagnostic flags in `benchmark_recombination.py` are all gated behind explicit booleans — leaving them on for the headline timing comparison contaminates the numbers (snapshot capture, etc., explicitly pause the wallclock, but instrument=True does add ~150 ns per moved mutation).
