#!/usr/bin/env python3
"""
Recombination Benchmarker

Benchmarks the native GRG recombination approach against the baseline
NumPy conversion approach (Recombine -> Dense Matrix).
"""

import os
import glob
import argparse
import time
import csv
import json
import gc
import sys
import platform
import tempfile
import numpy as np
from pathlib import Path
from dataclasses import dataclass, asdict
import pygrgl
import re
import psutil

# Import our custom modules (ensure these are in the same directory)
from grg_recombination import (
    simulate_grg_recombination,
    NonDuplicationRecombination,
    compute_grg_structural_stats,
)
from grg_numpy_baseline import grg_to_numpy, grg_to_numpy_parallel, estimate_numpy_memory, simulate_numpy_recombination
from multitree_check import compute_post_recomb_anc_counts, check_offspring, compute_liveness

# DEBUG Mode
debug = True


def get_process_memory_mb():
    """Returns the current process memory usage in MB."""
    process = psutil.Process(os.getpid())
    # rss = Resident Set Size (the actual RAM used)
    return process.memory_info().rss / (1024 * 1024)

# def get_breakpoints(bp_range, expected_crossovers=1.5):
#     num_bp = np.random.poisson(expected_crossovers)
#     if num_bp == 0:
#         return np.array([], dtype=int)

#     low, high = bp_range            # inclusive or exclusive? assuming inclusive high here
#     length = high - low + 1
#     num_bp = min(num_bp, length)

#     # Oversample a bit, then deduplicate
#     k = num_bp * 3                  # oversampling factor; small, since num_bp is tiny
#     candidates = np.random.randint(low, high + 1, size=k)
#     unique = np.unique(candidates)

#     if unique.size < num_bp:
#         # Very unlikely with small num_bp, but handle just in case
#         extra_needed = num_bp - unique.size
#         extra = np.random.randint(low, high + 1, size=extra_needed * 3)
#         unique = np.unique(np.concatenate([unique, extra]))

#     bp = np.sort(unique[:num_bp])
#     return bp


def recombination_intervals(h1, h2, bp, N):
    """
    Returns list of segments: [(source_parent_id, end_coord), ...]
    """
    # bp = self.get_breakpoints(N)
    # start = np.random.binomial(1, 0.5, 1)[0]
    start = h1
    parents = [h1, h2]
    
    segments = []
    for i, K in enumerate(bp):
        segments.append((parents[(start + i) % 2], K))
    segments.append((parents[(start + len(bp)) % 2], N))
    return segments

@dataclass
class BenchmarkResult:
    """Single benchmark measurement result"""
    file: str
    num_offspring: int
    implementation: str  # 'GRG Native' or 'NumPy Baseline'
    num_samples_initial: int
    num_snps: int
    num_runs: int
    num_bp: float
    mean_bp: float
    mean_time_ms: float
    std_time_ms: float
    min_time_ms: float
    max_time_ms: float
    nodes_added: int = 0
    edges_added: int = 0
    memory_mb: float = 0.0
    # Liveness diagnostic (averaged across timed_runs * num_generations snapshots).
    # Populated for GRG Native only; remains 0.0 for the NumPy baseline row.
    mean_dead_nodes: float = 0.0
    mean_dead_pct: float = 0.0
    mean_dead_mutations: float = 0.0
    mean_dead_edges: float = 0.0
    mean_dead_edges_pct: float = 0.0

@dataclass
class SystemInfo:
    """System and environment information"""
    numpy_version: str
    pygrgl_version: str
    platform: str
    cpu_count: int
    timestamp: str

class RecombinationBenchmarker:
    def __init__(
        self,
        grg_files_dir: Path,
        output_dir: Path,
        num_warmup: int = 1,
        num_runs: int = 3,
        num_offspring_per_couple: int = 2,
        num_generations: int = 2,
        memory_limit_mb: float = 15000.0,
        include_numpy: bool = True,
        verification: bool = False,
        run_diagnostics: bool = False,
        serialize: bool = False,
        profile: bool = False,
    ):
        self.grg_dir = grg_files_dir
        self.output_dir = output_dir
        self.num_warmup = num_warmup
        self.num_runs = num_runs
        self.num_offspring_per_couple = num_offspring_per_couple
        self.num_couples = None  # will be determined per file based on num_samples
        self.num_generations = num_generations
        self.memory_limit_mb = memory_limit_mb
        self.include_numpy = include_numpy
        self.verification = verification
        self.run_diagnostics = run_diagnostics
        self.serialize = serialize
        self.profile = profile
        self.results = []
        # Always-initialized so that verification_checks (which is a separate
        # concern from --diagnostics) can stash per-file results without
        # depending on whether structural-stats / diagnostic-pass ran.
        self.diagnostics = {}
        
        try:
            pygrgl_ver = pygrgl.__version__
        except AttributeError:
            pygrgl_ver = "unknown"
            
        import multiprocessing
        self.system_info = SystemInfo(
            numpy_version=np.__version__,
            pygrgl_version=pygrgl_ver,
            platform=f"{platform.system()}-{platform.release()}-{platform.machine()}",
            cpu_count=multiprocessing.cpu_count(),
            timestamp=time.strftime("%Y-%m-%d %H:%M:%S")
        )

    def _parse_expected_size(self, filename: str) -> tuple:
        """Extract expected individuals and SNPs from filename."""
        inds_match = re.search(r'(\d+)inds', filename)
        snps_match = re.search(r'(\d+)([kKmM])?_snps', filename)
        
        inds = int(inds_match.group(1)) if inds_match else 0
        
        snps = 0
        if snps_match:
            val = int(snps_match.group(1))
            mult = snps_match.group(2)
            if mult:
                mult = mult.lower()
                if mult == 'k': snps = val * 1000
                elif mult == 'm': snps = val * 1000000
            else:
                snps = val
                
        return inds, snps

    def print_header(self):
        print("=" * 80)
        print("GRG Recombination vs NumPy Baseline Benchmark")
        print("=" * 80)
        print(f"System: {self.system_info.platform}")
        print(f"CPU Cores: {self.system_info.cpu_count}")
        print(f"NumPy: {self.system_info.numpy_version}")
        print(f"Warmup runs: {self.num_warmup}, Timed runs: {self.num_runs}")
        print(f"Offspring per couple: {self.num_offspring_per_couple}")
        if self.include_numpy:
            print(
                "NumPy baseline: grg_numpy_baseline.simulate_numpy_recombination "
                "(one breakpoint draw per couple; sibling rows partition p1/p2 per segment)"
            )
        print(f"Memory limit: {self.memory_limit_mb} MB")
        print(f"Diagnostics: {'ON' if self.run_diagnostics else 'OFF'} "
              f"(structural stats + instrumented diagnostic pass + liveness/deadweight)")
        print(f"Verification: {'ON' if self.verification else 'OFF'} "
              f"(audit identities + multitree-violation check after each non-warmup run)")
        print(f"Save+Reload measurement: {'ON' if self.serialize else 'OFF'} "
              f"(per-generation save_grg+load_mutable_grg on final timed run, off the timed path; "
              f"passive -- graph is not replaced)")
        print(f"cProfile pass: {'ON' if self.profile else 'OFF'} "
              f"(one generation under cProfile, top functions by cumtime + tottime)")
        print("=" * 80)

    def get_breakpoints(self, bp_range, expected_crossovers=1.5):
        num_bp = np.random.poisson(expected_crossovers)
        if num_bp == 0:
            return np.array([], dtype=int), num_bp

        low, high = bp_range  # assuming exclusive low and high
        length = high - low + 1
        num_bp = min(num_bp, length)

        # Oversample a bit, then deduplicate
        k = num_bp * 3  # oversampling factor; small, since num_bp is tiny
        candidates = np.random.randint(low + 1, high, size=k)
        unique = np.unique(candidates)

        if unique.size < num_bp:
            # Very unlikely with small num_bp, but handle just in case
            extra_needed = num_bp - unique.size
            extra = np.random.randint(low + 1, high, size=extra_needed * 3)
            unique = np.unique(np.concatenate([unique, extra]))

        bp = np.sort(unique[:num_bp])
        return bp, num_bp

    def run_benchmarks(self):
        self.print_header()

        # 1. Use pathlib to find the files (returns real Path objects, not strings!)
        if self.grg_dir.is_dir():
            grg_files = list(self.grg_dir.glob("*.grg"))
        else:
            grg_files = [self.grg_dir] if self.grg_dir.suffix == ".grg" else []

        # 2. Sort from smallest to largest using Path's native stat().st_size
        grg_files.sort(key=lambda x: x.stat().st_size)
        
        if not grg_files:
            print(f"No .grg files found in {self.grg_dir}")
            return []

        print("Files to process (Smallest to Largest):")
        for f in grg_files:
            # f is now a real Path object, so .name works perfectly here!
            print(f"  {f.name}")
        print("=" * 80)

        # Set base csv output filename prefix (we'll append file-specific suffix if benchmark is run on 1 file)
        filename_prefix = "benchmark_recombination_results"
        if len(grg_files) == 1:
            # If only one file, include its name in the output filename prefix for clarity
            filename_prefix += f"_{grg_files[0].stem}"

        for file_path in grg_files:
            # file_path is a real Path object, so _benchmark_file will accept it perfectly!
            self._benchmark_file(file_path)
            
            self.save_results(self.output_dir, filename_prefix)
            print(f"[*] Checkpoint saved! Results safely written to disk.")
        return self.results

    def _benchmark_file(self, file_path: Path):
        print(f"\nProcessing: {file_path.name}")
        exp_inds, exp_snps = self._parse_expected_size(file_path.name)
        if exp_inds and exp_snps:
            print(f" Expected: {exp_inds} individuals, {exp_snps} SNPs")

        # Load once just to get base stats
        g_base = pygrgl.load_mutable_grg(str(file_path))
        base_samples = g_base.num_samples
        base_mutations = g_base.num_mutations
        base_nodes = g_base.num_nodes
        base_edges = g_base.num_edges
        # Shared physical interval for breakpoint draws (GRG + NumPy baseline).
        base_genome = g_base.bp_range
        base_couples = base_samples // 2
        self.num_couples = base_couples

        print(f" Genome (bp_range for recombination): {base_genome}")
        print(f" Couples: {base_couples} (from {base_samples} samples)")

        # ---- Structural stats (one-shot, before any recombination) ----
        # Gated by --diagnostics: compute & report structural fingerprint
        # (num_nodes, num_edges, roots, mutations-per-node distribution, etc.).
        # Off by default since this is for understanding GRG topology, not
        # for the timing/memory benchmark itself.
        if self.run_diagnostics:
            print("\nComputing structural stats...")
            struct_t0 = time.perf_counter()
            structural_stats = compute_grg_structural_stats(g_base)
            print(f"  done in {time.perf_counter() - struct_t0:.2f}s")
            for key in ("num_nodes", "num_edges", "num_samples", "num_mutations",
                        "num_roots", "internal_per_sample", "genome_length"):
                print(f"  {key:25s} {structural_stats[key]}")
            for key in ("mutations_per_node", "up_fanout", "sample_depth"):
                s = structural_stats[key]
                if s.get("count", 0):
                    print(f"  {key:25s} mean={s['mean']:.2f} p50={s['p50']} p95={s['p95']} max={s['max']}")
            self.diagnostics[file_path.name] = {"structural": structural_stats}

        # Estimate memory for the POST-recombination NumPy array
        # It will have (base_samples + num_offspring) columns
        est_mb = estimate_numpy_memory(base_mutations, base_samples + self.num_offspring_per_couple * base_couples)

        skip_numpy = not self.include_numpy
        if skip_numpy:
            print(f"  [Skipping NumPy] disabled via --skip-numpy")
        elif est_mb > self.memory_limit_mb:
            print(f"  [Skipping NumPy] Estimated memory ({est_mb:.1f} MB) exceeds limit ({self.memory_limit_mb} MB)")
            skip_numpy = True

        print(f" Actual Base: {base_samples} samples, {base_mutations} mutations, {base_nodes} nodes")
        print(f" Estimated memory for NumPy array (post-recomb): {est_mb:.1f} MB")

        # Free base graph
        del g_base
        gc.collect()

        # generate segments
        # run grg recombination
        # run numpy recombination
        # repeat steps 1-3 for num_runs and num_warmup

        if not skip_numpy:
            print(f"\nGenerating Base NumPy Population Matrix...")
            base_numpy_population_matrix = grg_to_numpy_parallel(pygrgl.load_mutable_grg(str(file_path)))

        grg_times = []
        grg_sizes = []
        grg_size_changes = []
        numpy_times = []
        total_nodes_added = 0
        total_edges_added = 0
        # Liveness diagnostic: one entry per (timed_run, generation). Off the
        # timed path -- we pause the wallclock timer around each call so the
        # reported GRG runtime excludes BFS-upward time.
        liveness_snapshots = []

        for i in range(self.num_warmup + self.num_runs):

            print(f"  [Run {i+1}] Loading graph from disk...")
            load_start = time.time()
            g = pygrgl.load_mutable_grg(str(file_path))
            print(f"  [Run {i+1}] Graph loaded in {time.time() - load_start:.2f} seconds.")

            # testcases = []

            # # Generate testcase breakpoints and segments for all offspring in this run (for both GRG and NumPy)
            # if debug:
            #     print(f"  [Run {i+1}] Generating testcases")
            # for gen in range(self.num_generations):
            #     gen_testcases = []
            #     for couple in range(base_couples):
            #         couple_testcases = []
            #         for offspring in range(4):  # Each pair produces 2 children and each child has 2 sets of breakpoints due to diploid instance
            #             print("testcase breakpoints:", couple, offspring)
            #             bp = get_breakpoints(base_genome)
            #             couple_testcases.append(bp)
            #         gen_testcases.append(couple_testcases)
            #     testcases.append(gen_testcases)
            # if debug:
            #     print(f"  [Run {i+1}] Testcases generated.")



            # GRG Native
            gc.collect()
            mem_before = get_process_memory_mb()
            print(f"\nStarting GRG Recombination Benchmarking...")
            if debug:
                print(f"\n  [Run {i+1}] Simulating {self.num_generations} generations with GRG recombination...")
            total_grg_bp = 0
            all_offspring_ids = []
            start = time.perf_counter()
            recomb = NonDuplicationRecombination(g)
            for gen in range(self.num_generations):
                # Run GRG recombination for this generation
                if debug:
                    print(f"\n    [Gen {gen+1}] Simulating generation {gen+1} with GRG recombination...")
                # pr = cProfile.Profile()
                # pr.enable()
                offspring_ids, gen_bp = simulate_grg_recombination(self, recomb, base_genome, N=base_genome[1])
                # pr.disable()
                # pstats.Stats(pr).sort_stats('cumtime').print_stats(25)  # Print top 25 cumulative time functions
                all_offspring_ids.extend(offspring_ids)
                if i >= self.num_warmup:
                    total_grg_bp += gen_bp
                    # Liveness snapshot off the timed path. Gated on --diagnostics
                    # so production benchmarks don't pay even the (cheap) O(|V|+|E|)
                    # BFS-upward cost when the user isn't asking for it.
                    # simulate_grg_recombination calls set_samples(new_offspring_ids)
                    # at its tail, so the GRG's current sample set IS the new
                    # generation; compute_liveness uses those as BFS-up seeds. We
                    # pause the wallclock around the call so `elapsed` excludes the
                    # BFS time even when the flag is on.
                    if self.run_diagnostics:
                        t_pause = time.perf_counter()
                        snap = compute_liveness(g)
                        snap['run_index'] = i
                        snap['gen'] = gen
                        liveness_snapshots.append(snap)
                        start += (time.perf_counter() - t_pause)

                # ----- Per-generation save_grg + load_mutable_grg measurement -----
                # Off the timed path; pauses the wallclock around the call so
                # `elapsed` excludes save+load time. Gated on (last timed iteration)
                # so per-file we get one measurement per generation, all from the
                # same accumulating un-trimmed graph state.
                # NOTE: this is a passive measurement -- g is NOT replaced with the
                # reloaded (trimmed) version. Subsequent generations continue on
                # the original un-trimmed graph. To actually run subsequent gens
                # against the trimmed graph would require rebuilding recomb's
                # caches against the renumbered post-reload node IDs.
                if (self.serialize
                        and i == (self.num_warmup + self.num_runs - 1)):
                    t_pause = time.perf_counter()
                    print(f"\n    [Gen {gen+1}] Running save_grg + load_mutable_grg...")
                    pre_nodes = g.num_nodes
                    pre_edges = g.num_edges
                    with tempfile.NamedTemporaryFile(suffix='.grg', delete=False) as _tf:
                        sr_tmp_path = _tf.name
                    try:
                        t0 = time.perf_counter()
                        pygrgl.save_grg(g, sr_tmp_path)
                        save_s = time.perf_counter() - t0
                        file_size_mb = os.path.getsize(sr_tmp_path) / (1024 * 1024)

                        t0 = time.perf_counter()
                        g_reloaded = pygrgl.load_mutable_grg(sr_tmp_path)
                        load_s = time.perf_counter() - t0

                        post_nodes = g_reloaded.num_nodes
                        post_edges = g_reloaded.num_edges
                        node_delta = pre_nodes - post_nodes
                        edge_delta = pre_edges - post_edges
                        node_pct = 100.0 * node_delta / max(1, pre_nodes)
                        edge_pct = 100.0 * edge_delta / max(1, pre_edges)
                        total_s = save_s + load_s

                        print(f"    [Save+Reload gen {gen+1}] pre:   {pre_nodes:>12,} nodes  {pre_edges:>12,} edges")
                        print(f"    [Save+Reload gen {gen+1}] post:  {post_nodes:>12,} nodes  {post_edges:>12,} edges")
                        print(f"    [Save+Reload gen {gen+1}] delta: -{node_delta:>11,} nodes (-{node_pct:5.2f}%)  "
                              f"-{edge_delta:>11,} edges (-{edge_pct:5.2f}%)")
                        print(f"    [Save+Reload gen {gen+1}] save: {save_s:.2f}s  load: {load_s:.2f}s  "
                              f"total: {total_s:.2f}s  file: {file_size_mb:.1f} MB")

                        self.diagnostics.setdefault(file_path.name, {}) \
                            .setdefault('save_reload_per_gen', []).append({
                                'gen': gen,
                                'pre_nodes': pre_nodes,
                                'pre_edges': pre_edges,
                                'post_nodes': post_nodes,
                                'post_edges': post_edges,
                                'node_reduction_pct': node_pct,
                                'edge_reduction_pct': edge_pct,
                                'save_s': save_s,
                                'load_s': load_s,
                                'total_s': total_s,
                                'file_size_mb': file_size_mb,
                                'measured_on_iteration': i,
                            })
                        del g_reloaded
                    finally:
                        if os.path.exists(sr_tmp_path):
                            os.remove(sr_tmp_path)
                    gc.collect()
                    start += time.perf_counter() - t_pause

                if debug:
                    print(f"    [Gen {gen+1}] Generation's Breakpoints: {gen_bp}. Total so far: {total_grg_bp}")
            elapsed = time.perf_counter() - start
            if debug:
                print(f"\n  [Run {i+1}] Simulation finished in {elapsed:.4f} seconds.")

            if i >= self.num_warmup:
                grg_times.append(elapsed * 1000)  # Convert to ms

            # Grab node and edge counts on timed runs for space complexity metrics
            if i >= self.num_warmup:
                total_nodes_added += g.num_nodes - base_nodes
                total_edges_added += g.num_edges - base_edges

            # ----- Verification block -----
            # Off the timed path; runs only on non-warmup runs when opted in
            # via --verification. Two layers:
            #   (1) Three audit-identity checks (cheap; just summing audit
            #       counters that recomb tracked during the run).
            #   (2) Multitree-violation check (Approach B; per-offspring BFS,
            #       expensive on biobank-scale data).
            if self.verification and i >= self.num_warmup:
                print(f"\n  [Run {i+1}] Running verification checks...")

                # (1) Audit identities. recomb.audit accumulates across all
                # generations of this run (we never reset_stats between gens
                # on the non-instrumented path).
                audit_results = recomb.audit_check(raise_on_fail=False)
                audit_pass = sum(1 for r in audit_results.values() if r['pass'])
                print(f"  [Run {i+1}] Audit identities: {audit_pass}/"
                      f"{len(audit_results)} pass")
                for r in audit_results.values():
                    mark = "OK" if r['pass'] else "FAIL"
                    print(f"    [{mark}] {r['desc']}: lhs={r['lhs']} rhs={r['rhs']}")

                # (2) Multitree cardinality check.
                mt_start = time.perf_counter()
                anc_count = compute_post_recomb_anc_counts(g)
                mt_result = check_offspring(g, all_offspring_ids, anc_count)
                mt_elapsed = time.perf_counter() - mt_start
                mt_sum = mt_result['summary']
                print(f"  [Run {i+1}] Multitree check done in {mt_elapsed:.2f}s: "
                      f"{mt_sum['violating']:,} / {mt_sum['total_checked']:,} "
                      f"offspring violate ({100*mt_sum['violation_rate']:.2f}%); "
                      f"{mt_sum['total_redundant_paths']:,} redundant paths total")

                self.diagnostics.setdefault(file_path.name, {}) \
                    .setdefault('verification_checks', []).append({
                        'run_index': i,
                        'is_warmup': False,
                        'audit_identities': audit_results,
                        'multitree': {
                            'wallclock_s': mt_elapsed,
                            'offspring_checked_across_gens': len(all_offspring_ids),
                            'summary': mt_sum,
                            'first_example': mt_result['first_example'],
                        },
                    })

            gc.collect()
            mem_after = get_process_memory_mb()
            grg_sizes.append(mem_after)
            grg_size_changes.append(mem_after - mem_before)
            print(f"  [Run {i+1}] Memory usage: {mem_after:.1f} MB (Delta: +{mem_after - mem_before:.1f} MB)")

            del g
            del recomb

            gc.collect()
            # Numpy
            if not skip_numpy:
                print(f"\nStarting NumPy Recombination Benchmarking...")
                if debug:
                    print(f"\n  [Run {i+1}] Simulating {self.num_generations} generations with NumPy recombination...")
                offspring_matrix = base_numpy_population_matrix.copy()
                total_numpy_bp = 0
                start = time.perf_counter()
                for gen in range(self.num_generations):
                    # Run NumPy recombination for this generation
                    if debug:
                        print(f"\n    [Gen {gen+1}] Simulating generation {gen+1} with NumPy recombination...")
                    # self provides get_breakpoints(bp_range, ...) and num_offspring_per_couple
                    # (see grg_numpy_baseline.simulate_numpy_recombination).
                    offspring_matrix, gen_bp = simulate_numpy_recombination(
                        self,
                        offspring_matrix,
                        bp_range=base_genome,
                        expected_crossovers=1.5,
                    )
                    if i >= self.num_warmup:
                        total_numpy_bp += gen_bp
                    if debug:
                        print(f"    [Gen {gen+1}] Generation's Breakpoints: {gen_bp}. Total so far: {total_numpy_bp}")
                elapsed = time.perf_counter() - start
                if debug:
                    print(f"\n  [Run {i+1}] Simulation finished in {elapsed:.4f} seconds.")
                if i >= self.num_warmup:
                    numpy_times.append(elapsed * 1000)  # Convert to ms

        # ----- Diagnostic pass with full instrumentation -----
        # Gated by --diagnostics: one additional GRG run with `instrument=True`
        # so we capture phase breakdowns, C++ call costs (esp. add_mutation /
        # remove_mutation), per-offspring + per-bubble distributions per
        # generation, and the full audit-1 decision-case histogram.
        # Wallclock here is NOT used for the headline numbers above.
        if self.run_diagnostics:
            print(f"\nRunning instrumented diagnostic pass...")
            diag_g = pygrgl.load_mutable_grg(str(file_path))
            diag_recomb = NonDuplicationRecombination(diag_g, instrument=True)
            # _build_ancestral_caches runs in __init__ and writes
            # init_caches_time into stats. Capture it now -- the first
            # reset_stats() below would zero it out, and the init pass is
            # a one-time construction cost that doesn't belong in per-gen
            # snapshots anyway.
            init_caches_time_s = diag_recomb.stats["init_caches_time"]
            print(f"  [Diag] One-pass cache init: {init_caches_time_s:.2f}s "
                  f"(amortized across {self.num_generations} gen(s))")
            per_generation_stats = []
            diag_total_start = time.perf_counter()
            for gen in range(self.num_generations):
                print(f"  [Diag] Generation {gen+1}/{self.num_generations}...")
                diag_recomb.reset_stats()  # also resets audit
                gen_start = time.perf_counter()
                simulate_grg_recombination(self, diag_recomb, base_genome, N=base_genome[1])
                gen_elapsed = time.perf_counter() - gen_start
                snapshot = dict(diag_recomb.stats)
                # ---- Audit 1: snapshot per-generation case histogram ----
                # `dict(...)` makes a shallow copy so subsequent reset_stats()
                # doesn't clobber the saved values.
                snapshot["audit"] = dict(diag_recomb.audit)
                snapshot["generation_index"] = gen
                snapshot["generation_wallclock_s"] = gen_elapsed
                snapshot["num_nodes_after_gen"] = diag_g.num_nodes
                snapshot["num_edges_after_gen"] = diag_g.num_edges
                # Total edges added this generation = all connect() calls. The audit
                # tracks both direct-attach and extract-bubble connects;
                # stats['connect_calls'] only covers extract-bubble (timed path).
                edges_added_this_gen = (
                    snapshot["audit"]["connect_calls_in_attach"]
                    + snapshot["audit"]["connect_calls_in_extract"]
                )
                snapshot["edges_added_this_gen"] = edges_added_this_gen
                per_generation_stats.append(snapshot)
                print(f"    wallclock={gen_elapsed:.2f}s "
                      f"offspring={snapshot['offspring_count']} "
                      f"bubbles={snapshot['bubbles_created']} "
                      f"muts_moved={snapshot['mutations_moved']} "
                      f"visits={snapshot['visits_total']} "
                      f"edges_added={edges_added_this_gen}")
                if snapshot["remove_mutation_calls"]:
                    rm_us = snapshot["remove_mutation_time"] / snapshot["remove_mutation_calls"] * 1e6
                    print(f"    remove_mutation per-call: {rm_us:.2f} us "
                          f"(total {snapshot['remove_mutation_time']:.2f}s "
                          f"over {snapshot['remove_mutation_calls']:,} calls)")
                if snapshot["add_mutation_calls"]:
                    am_us = snapshot["add_mutation_time"] / snapshot["add_mutation_calls"] * 1e6
                    print(f"    add_mutation per-call:    {am_us:.2f} us")
            diag_total = time.perf_counter() - diag_total_start
            print(f"  Diagnostic pass total: {diag_total:.2f}s")

            # ---- Audit 1: aggregate per-gen audits and print histogram ----
            # Each per-gen snapshot has its own audit dict (because reset_stats
            # zeroes them between gens). Sum across gens to get the full-run
            # histogram, run audit_check on it, and pretty-print.
            aggregated_audit = NonDuplicationRecombination._fresh_audit()
            for snap in per_generation_stats:
                for k, v in snap.get("audit", {}).items():
                    aggregated_audit[k] = aggregated_audit.get(k, 0) + v
            diag_recomb.audit_summary(
                audit=aggregated_audit,
                header=f"{file_path.name} -- {self.num_generations} gen(s) cumulative",
            )

            del diag_g, diag_recomb
            gc.collect()

            # Stash per-generation snapshots on the diagnostic record for this file.
            # `aggregated_audit` is also stashed so consumers of the JSON don't
            # have to re-aggregate from per_generation entries.
            self.diagnostics.setdefault(file_path.name, {})["per_generation"] = per_generation_stats
            self.diagnostics[file_path.name]["audit_aggregated"] = aggregated_audit
            self.diagnostics[file_path.name]["diagnostic_total_wallclock_s"] = diag_total
            # One-time __init__ cost for the topological cache prebuild
            # (separate from per_generation since it doesn't recur per gen).
            self.diagnostics[file_path.name]["init_caches_time_s"] = init_caches_time_s

        # ----- cProfile pass (one generation, fresh graph) -----
        # Gated by --profile. Uses cProfile rather than manual perf_counter probes
        # inside the hot loop, so there is no per-iteration overhead. The profiler
        # accounts for time at function-call granularity, giving a breakdown across:
        #   _recurse_attach (tottime = pure traversal; cumtime includes callees)
        #   _get_node_and_ancestor_span  (ancestral interval DFS on cache miss)
        #   _get_ancestral_coverage      (anc-cov cache miss path)
        #   _get_up_edges_cached         (edge fetching)
        #   _get_node_mutations          (mutation cache miss path)
        #   _extract_bubble              (bubble creation)
        #   _apply_pending_bubbles       (deferred mutation moves)
        #   grg.connect / grg.make_node  (C++ boundary calls)
        # Wall-clock here is NOT used for headline numbers.
        if self.profile:
            import cProfile
            import pstats
            import io

            print(f"\nRunning cProfile pass (1 generation, fresh graph)...")
            prof_g = pygrgl.load_mutable_grg(str(file_path))
            prof_recomb = NonDuplicationRecombination(prof_g)

            pr = cProfile.Profile()
            pr.enable()
            simulate_grg_recombination(self, prof_recomb, base_genome, N=base_genome[1])
            pr.disable()

            stream = io.StringIO()
            ps = pstats.Stats(pr, stream=stream).strip_dirs()

            ps.sort_stats('cumtime').print_stats(30)
            cumtime_output = stream.getvalue()
            stream.seek(0); stream.truncate(0)

            ps.sort_stats('tottime').print_stats(30)
            tottime_output = stream.getvalue()

            print(f"\n--- cProfile: top 30 by cumulative time ---")
            print(cumtime_output)
            print(f"--- cProfile: top 30 by total time (self only, excludes callees) ---")
            print(tottime_output)

            self.diagnostics.setdefault(file_path.name, {})['profile'] = {
                'cumtime_top30': cumtime_output,
                'tottime_top30': tottime_output,
            }

            del prof_g, prof_recomb
            gc.collect()

        # ----- Liveness aggregation across (timed_runs * num_generations) snapshots -----
        # Per-snapshot keys we average. dead_samples is a sanity counter (always 0 if
        # the algorithm is correct); included so a non-zero average loudly fails.
        # Aggregation + diagnostics stash run only when --diagnostics is on; when off,
        # liveness_means stays zero-valued so the BenchmarkResult plumb below records
        # 0.0 for the mean_dead_* columns (CSV stays well-formed).
        liveness_keys = (
            'total_nodes', 'alive', 'dead', 'dead_pct',
            'num_samples', 'dead_samples', 'dead_roots',
            'dead_with_mutations', 'dead_internal_empty', 'dead_mutation_count',
            'total_edges', 'alive_alive_edges', 'dead_alive_edges',
            'dead_dead_edges', 'dead_edges', 'dead_edges_pct',
        )
        liveness_means = {k: 0.0 for k in liveness_keys}
        liveness_by_gen = []

        if self.run_diagnostics and liveness_snapshots:
            n_snaps = len(liveness_snapshots)
            liveness_means = {
                k: sum(s[k] for s in liveness_snapshots) / n_snaps
                for k in liveness_keys
            }
            # Per-generation averages (across timed runs, for each generation index).
            by_gen = {}
            for s in liveness_snapshots:
                by_gen.setdefault(s['gen'], []).append(s)
            for g_idx in sorted(by_gen):
                snaps = by_gen[g_idx]
                liveness_by_gen.append({
                    'gen': g_idx,
                    'snapshot_count': len(snaps),
                    **{f'mean_{k}': sum(s[k] for s in snaps) / len(snaps) for k in liveness_keys},
                })

            self.diagnostics.setdefault(file_path.name, {})['liveness'] = {
                'snapshot_count': n_snaps,
                'snapshots': liveness_snapshots,
                'means_over_all_snapshots': liveness_means,
                'means_by_generation': liveness_by_gen,
            }

        grg_mean, grg_std, grg_sizes_mean, grg_size_changes_mean = np.mean(grg_times), np.std(grg_times), np.mean(grg_sizes), np.mean(grg_size_changes)
        grg_bp_mean = total_grg_bp / (self.num_runs * self.num_generations)
        if not skip_numpy:
            numpy_mean, numpy_std = np.mean(numpy_times), np.std(numpy_times)
            numpy_bp_mean = total_numpy_bp / (self.num_runs * self.num_generations)
        actual_offspring_generated = self.num_couples * self.num_offspring_per_couple * self.num_generations
        nodes_added_per_run = total_nodes_added / self.num_runs
        edges_added_per_run = total_edges_added / self.num_runs

        self.results.append(BenchmarkResult(
            file=file_path.name, num_offspring=actual_offspring_generated,
            implementation="GRG Native", num_samples_initial=base_samples,
            num_snps=base_mutations, num_runs=self.num_runs,
            num_bp=total_grg_bp, mean_bp=grg_bp_mean,
            mean_time_ms=grg_mean, std_time_ms=grg_std,
            min_time_ms=np.min(grg_times), max_time_ms=np.max(grg_times),
            nodes_added=nodes_added_per_run,
            edges_added=edges_added_per_run,
            memory_mb=grg_sizes_mean,
            mean_dead_nodes=liveness_means['dead'],
            mean_dead_pct=liveness_means['dead_pct'],
            mean_dead_mutations=liveness_means['dead_mutation_count'],
            mean_dead_edges=liveness_means['dead_edges'],
            mean_dead_edges_pct=liveness_means['dead_edges_pct'],
        ))

        if not skip_numpy:
            self.results.append(BenchmarkResult(
                file=file_path.name, num_offspring= actual_offspring_generated,
                implementation="NumPy Baseline", num_samples_initial=base_samples,
                num_snps=base_mutations, num_runs=self.num_runs,
                num_bp=total_numpy_bp, mean_bp=numpy_bp_mean,
                mean_time_ms=numpy_mean, std_time_ms=numpy_std,
                min_time_ms=np.min(numpy_times), max_time_ms=np.max(numpy_times),
                nodes_added=0, edges_added=0, memory_mb=est_mb
            ))

        print(f"\nResults for {file_path.name}:\n")

        print(f"  GRG Native:     {grg_mean:.2f}ms ± {grg_std:.2f}ms")
        print(f"  Space Delta:    +{nodes_added_per_run:.2f} nodes created per run on average")
        print(f"  Edge Delta:     +{edges_added_per_run:.2f} edges created per run on average")
        print(f"  GRG Breakpoints: {grg_bp_mean:.2f} average breakpoints per generation")
        print(f"  GRG Memory:     {grg_sizes_mean:.1f} MB average resident memory")
        print(f"  Memory Delta:   +{grg_size_changes_mean:.1f} MB (GRG recombination memory increase), ")

        # ----- Liveness (deadweight) averages -----
        # Gated on --diagnostics; the snapshot-capture loop above no-ops when off,
        # so liveness_snapshots is empty and we silently skip the print block.
        if self.run_diagnostics and liveness_snapshots:
            n_snaps = len(liveness_snapshots)
            print(f"\n  Liveness (avg over {self.num_runs} timed runs * "
                  f"{self.num_generations} gens = {n_snaps} snapshots):")
            print(f"    mean total nodes (post-gen):  {liveness_means['total_nodes']:.2f}")
            print(f"    mean alive nodes:             {liveness_means['alive']:.2f}")
            print(f"    mean dead nodes:              {liveness_means['dead']:.2f}  "
                  f"({liveness_means['dead_pct']:.1f}%)")
            print(f"    mean dead w/ mutations:       {liveness_means['dead_with_mutations']:.2f}")
            print(f"    mean dead empty internals:    {liveness_means['dead_internal_empty']:.2f}")
            print(f"    mean dead roots:              {liveness_means['dead_roots']:.2f}")
            print(f"    mean orphaned mutations:      {liveness_means['dead_mutation_count']:.2f}")
            print(f"    mean total edges (post-gen):  {liveness_means['total_edges']:.2f}")
            print(f"    mean dead edges:              {liveness_means['dead_edges']:.2f}  "
                  f"({liveness_means['dead_edges_pct']:.1f}%)")
            print(f"      dead-dead (both ends dead): {liveness_means['dead_dead_edges']:.2f}")
            print(f"      dead-alive (dangling):      {liveness_means['dead_alive_edges']:.2f}")
            if liveness_means['dead_samples']:
                print(f"    [!] mean dead_samples={liveness_means['dead_samples']:.2f} "
                      f"(should be 0 -- algorithm sanity violation)")
            if liveness_by_gen:
                print(f"\n    Per-generation averages across {self.num_runs} timed runs:")
                for entry in liveness_by_gen:
                    print(f"      gen {entry['gen']+1}: "
                          f"dead={entry['mean_dead']:.2f}  "
                          f"pct={entry['mean_dead_pct']:.1f}%  "
                          f"orphan_muts={entry['mean_dead_mutation_count']:.2f}")

        # ----- Save+Reload per-generation recap -----
        # Brief restatement of each per-generation save_grg+load_mutable_grg
        # measurement (full numbers printed inline during the last timed run).
        sr_list = self.diagnostics.get(file_path.name, {}).get('save_reload_per_gen')
        if sr_list:
            print(f"\n  Save+Reload (per-generation, last timed run, "
                  f"{len(sr_list)} measurement(s)):")
            for sr in sr_list:
                print(f"    gen {sr['gen']+1}: "
                      f"nodes {sr['pre_nodes']:>12,} -> {sr['post_nodes']:>12,} "
                      f"(-{sr['node_reduction_pct']:5.2f}%)  "
                      f"edges {sr['pre_edges']:>12,} -> {sr['post_edges']:>12,} "
                      f"(-{sr['edge_reduction_pct']:5.2f}%)  "
                      f"cost {sr['total_s']:.2f}s  ({sr['file_size_mb']:.1f} MB)")
        print()

        if not skip_numpy:
            print(f"  NumPy Baseline: {numpy_mean:.2f}ms ± {numpy_std:.2f}ms")
            print(
                f"  NumPy Breakpoints: {numpy_bp_mean:.2f} mean total breakpoint count per generation "
                f"(one draw per mating couple)"
            )
            if grg_mean > 0:
                speedup = numpy_mean / grg_mean
                print(f"  Speedup:        {speedup:.2f}x (NumPy / GRG)")
            print(f"  Memory Footprint: {est_mb:.1f} MB dense matrix allocated")
        
        # # ---------------------------------------------------------
        # # Phase A: GRG Native Recombination
        # # ---------------------------------------------------------

        # print("Starting GRG Benchmarking...")

        # grg_times = []
        # nodes_added = 0
        
        # for i in range(self.num_warmup + self.num_runs):
        #     # We MUST reload the graph fresh every run so it doesn't infinitely grow
        #     print(f"  [Run {i+1}] Loading graph from disk...")
        #     load_start = time.time()
        #     g = pygrgl.load_mutable_grg(str(file_path))
        #     print(f"  [Run {i+1}] Graph loaded in {time.time() - load_start:.2f} seconds.")
            
        #     print(f"  [Run {i+1}] Simulating {self.num_offspring} offspring...")
        #     start = time.perf_counter()
        #     simulate_generation(g)
        #     elapsed = time.perf_counter() - start
        #     print(f"  [Run {i+1}] Simulation finished in {elapsed:.4f} seconds.")
            
        #     if i >= self.num_warmup:
        #         grg_times.append(elapsed * 1000)  # Convert to ms
                
        #     # Grab node count on the last run for space complexity metrics
        #     if i == (self.num_warmup + self.num_runs - 1):
        #         nodes_added = g.num_nodes - base_nodes
                
        #     del g
        #     gc.collect()

        # grg_mean, grg_std = np.mean(grg_times), np.std(grg_times)

        # actual_offspring_generated = (g_base.num_samples // 2) * 2
        
        # self.results.append(BenchmarkResult(
        #     file=file_path.name, num_offspring=actual_offspring_generated,
        #     implementation="GRG Native", num_samples_initial=base_samples,
        #     num_snps=base_mutations, num_runs=self.num_runs,
        #     mean_time_ms=grg_mean, std_time_ms=grg_std,
        #     min_time_ms=np.min(grg_times), max_time_ms=np.max(grg_times),
        #     nodes_added=nodes_added, memory_mb=0.0
        # ))
        
        # print(f"  GRG Native:     {grg_mean:.2f}ms ± {grg_std:.2f}ms")
        # print(f"  Space Delta:    +{nodes_added} nodes created")

        # # ---------------------------------------------------------
        # # Phase B: NumPy Baseline
        # # ---------------------------------------------------------

        # print("Starting NumPy Baseline Benchmarking...")

        # # 1. SETUP (Not timed) - Convert GRG to Dense Matrix once
        # # This creates the starting population matrix (num_samples x num_snps)
        # numpy_population_matrix = grg_to_numpy(g)
        
        # numpy_times = []
        
        # for run in range(self.num_runs):
        #     gc.collect()
        #     start_time = time.perf_counter()
            
        #     # 2. EXECUTION (Timed) - Perform array-based recombination
        #     offspring_matrix = simulate_generation_numpy(numpy_population_matrix)
            
        #     end_time = time.perf_counter()
        #     numpy_times.append((end_time - start_time) * 1000)

        # np_mean = np.mean(numpy_times)
        # print(f"  NumPy Baseline: {np_mean:.2f}ms ± {np.std(numpy_times):.2f}ms")
        
        # # Calculate total memory of the baseline (initial matrix + offspring matrix)
        # total_np_memory = (numpy_population_matrix.nbytes + offspring_matrix.nbytes) / (1024 * 1024)
        # print(f"  Memory Footprint: {total_np_memory:.2f} MB dense matrices allocated")

        # if not self.include_numpy:
        #     return

        # if est_mb > self.memory_limit_mb:
        #     print(f"  [Skipping NumPy] Estimated memory ({est_mb:.1f} MB) exceeds limit ({self.memory_limit_mb} MB)")
        #     return

        # np_times = []
        # for i in range(self.num_warmup + self.num_runs):
        #     g = pygrgl.load_mutable_grg(str(file_path))
            
        #     # The entire baseline process: Recombine + Materialize to Dense Matrix
        #     start = time.perf_counter()
        #     simulate_generation(g, self.num_offspring)
        #     _ = grg_to_numpy(g)
        #     elapsed = time.perf_counter() - start
            
        #     if i >= self.num_warmup:
        #         np_times.append(elapsed * 1000)
                
        #     del g
        #     gc.collect()

        # np_mean, np_std = np.mean(np_times), np.std(np_times)
        
        # self.results.append(BenchmarkResult(
        #     file=file_path.name, num_offspring=self.num_offspring,
        #     implementation="NumPy Baseline", num_samples_initial=base_samples,
        #     num_snps=base_mutations, num_runs=self.num_runs,
        #     mean_time_ms=np_mean, std_time_ms=np_std,
        #     min_time_ms=np.min(np_times), max_time_ms=np.max(np_times),
        #     nodes_added=0, memory_mb=est_mb
        # ))

        # print(f"  NumPy Baseline: {np_mean:.2f}ms ± {np_std:.2f}ms")
        # if grg_mean > 0:
        #     speedup = np_mean / grg_mean
        #     print(f"  Speedup:        {speedup:.2f}x (NumPy / GRG)")
        # print(f"  Memory Footprint: {est_mb:.1f} MB dense matrix allocated")


    def save_results(self, output_dir: Path, filename_prefix: str = "benchmark_recombination_results"):
        output_dir.mkdir(parents=True, exist_ok=True)
        
        csv_path = output_dir / f"{filename_prefix}.csv"
        with open(csv_path, 'w', newline='') as f:
            if not self.results: return
            writer = csv.DictWriter(f, fieldnames=asdict(self.results[0]).keys())
            writer.writeheader()
            for res in self.results:
                writer.writerow(asdict(res))
                
        json_path = output_dir / f"{filename_prefix}.json"
        with open(json_path, 'w') as f:
            json.dump({
                "system_info": asdict(self.system_info),
                "config": {
                    "warmup_runs": self.num_warmup,
                    "timed_runs": self.num_runs,
                    "num_offspring": self.num_offspring_per_couple * self.num_couples * self.num_generations,
                    "memory_limit_mb": self.memory_limit_mb,
                    "num_generations": self.num_generations,
                    "num_offspring_per_couple": self.num_offspring_per_couple,
                    "profile": self.profile,
                },
                "results": [asdict(r) for r in self.results],
                "diagnostics": getattr(self, "diagnostics", {}),
            }, f, indent=2)

        print(f"\nResults saved to {csv_path} and {json_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GRG Recombination Benchmarker")
    parser.add_argument("--grg-dir", type=Path, default=Path("./grg_files"), help="Directory with GRG files")
    parser.add_argument("--output-dir", type=Path, default=Path("."), help="Directory to save benchmark results")
    parser.add_argument("--warmup", type=int, default=1, help="Warmup runs (default: 1)")
    parser.add_argument("--runs", type=int, default=3, help="Timed runs (default: 3)")
    parser.add_argument(
        "--offspring-per-couple",
        type=int,
        default=2,
        help="Offspring rows produced per mating couple each generation (NumPy: one shared breakpoint draw per couple; default: 2)",
    )
    parser.add_argument("--num-generations", type=int, default=2, help="Number of sequential recombination generations to simulate (default: 2)")
    parser.add_argument("--memory-limit", type=float, default=15000.0, help="Memory limit in MB (default: 15000)")
    parser.add_argument("--skip-numpy", action="store_true",
                        help="Skip the NumPy baseline entirely; no dense matrix is built")
    parser.add_argument("--verification", action="store_true",
                        help="After each non-warmup GRG run, run two verification layers: "
                             "(1) the three audit-identity checks (bubble_identity, "
                             "connect_identity, make_node_identity) on the live recomb's "
                             "audit counters; and (2) the multitree-violation check "
                             "(Approach B cardinality test). Off the timed path but adds "
                             "significant wallclock per file (e.g. ~2 min on 7500inds_1M). "
                             "Default off.")
    parser.add_argument("--diagnostics", action="store_true",
                        help="Run GRG-structural diagnostics for understanding the graph "
                             "(orthogonal to the timing/memory benchmark itself): pre-recomb "
                             "structural fingerprint (num_roots, mutations-per-node, fanout, "
                             "sample_depth distributions) and an extra fully-instrumented "
                             "diagnostic recomb pass (phase breakdowns, per-call C++ costs, "
                             "audit-1 decision-case histogram). Default off.")
    parser.add_argument("--serialize", action="store_true",
                        help="On the last timed run, after each generation, do a "
                             "pygrgl.save_grg(g, tmp) + pygrgl.load_mutable_grg(tmp) and "
                             "report pre/post (nodes, edges) plus save+load wallclock. "
                             "Off the timed path; runs num_generations times per file. The "
                             "graph itself is NOT replaced with the reloaded version -- "
                             "subsequent generations continue on the un-trimmed graph. Use "
                             "this to gauge per-generation how much the save-time simplifier "
                             "(extraneous-node strip + CSR repack) shrinks the accumulating "
                             "post-recombination graph. Default off.")
    parser.add_argument("--profile", action="store_true",
                        help="After timed runs, load the graph fresh and run one generation "
                             "under cProfile. Prints top 30 functions by cumulative time "
                             "(cumtime = function + callees) and by self time (tottime = "
                             "function body only, excluding callees). cumtime identifies the "
                             "costliest call chains; tottime isolates where CPU is actually "
                             "spent (e.g. pure traversal in _recurse_attach vs cache-miss DFS "
                             "in _get_node_and_ancestor_span). Off the timed path. Default off.")

    args = parser.parse_args()

    benchmarker = RecombinationBenchmarker(
        grg_files_dir=args.grg_dir,
        output_dir=args.output_dir,
        num_warmup=args.warmup,
        num_runs=args.runs,
        num_offspring_per_couple=args.offspring_per_couple,
        num_generations=args.num_generations,
        memory_limit_mb=args.memory_limit,
        include_numpy=not args.skip_numpy,
        verification=args.verification,
        run_diagnostics=args.diagnostics,
        serialize=args.serialize,
        profile=args.profile,
    )
    
    benchmarker.run_benchmarks()