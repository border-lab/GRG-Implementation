#!/usr/bin/env python3
"""
Single-run benchmark executor.

Executes exactly one run of one method (GRG or NumPy) on one .grg file and
writes a single-run result JSON.  Designed to be called many times in parallel
(e.g. via SLURM array jobs) with results aggregated afterwards by
benchmark_aggregate.py.

Example
-------
# GRG run:
python benchmark_run.py \\
    --grg-file ../grg_files/4000inds_500k_snps.grg \\
    --method grg \\
    --run-index 0 \\
    --num-generations 1 \\
    --offspring-per-couple 2 \\
    --output-dir ./output/runs

# NumPy run:
python benchmark_run.py \\
    --grg-file ../grg_files/4000inds_500k_snps.grg \\
    --method numpy \\
    --run-index 0 \\
    --num-generations 1 \\
    --offspring-per-couple 2 \\
    --output-dir ./output/runs
"""

import os
import sys
import gc
import json
import time
import argparse
import platform
import tempfile
import re
import threading
import numpy as np
from pathlib import Path

import pygrgl
import psutil

from grg_recombination import simulate_grg_recombination, compute_grg_structural_stats

if os.environ.get("GRG_BACKEND", "python").lower() in ("cpp", "c++", "native"):
    from grg_recombination_native import NonDuplicationRecombination
    _BACKEND = "cpp"
else:
    from grg_recombination import NonDuplicationRecombination
    _BACKEND = "python"

from grg_numpy_baseline import (
    grg_to_numpy_parallel,
    estimate_numpy_memory,
    simulate_numpy_recombination,
)
from multitree_check import compute_post_recomb_anc_counts, check_offspring, compute_liveness
from haplotype_oracle import check_offspring_haplotypes


def get_process_memory_mb():
    return psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024)


class PeakRSSTracker:
    """Background thread that polls process RSS to capture peak memory usage.

    Usage:
        tracker = PeakRSSTracker(interval=0.1)
        tracker.start()
        # ... work ...
        tracker.stop()
        print(tracker.peak_mb, tracker.baseline_mb)
    """

    def __init__(self, interval=0.1):
        self._interval = interval
        self._process = psutil.Process(os.getpid())
        self._peak = 0.0
        self._baseline = 0.0
        self._stop = threading.Event()
        self._thread = None

    def start(self):
        self._baseline = self._process.memory_info().rss / (1024 * 1024)
        self._peak = self._baseline
        self._stop.clear()
        self._thread = threading.Thread(target=self._poll, daemon=True)
        self._thread.start()

    def _poll(self):
        while not self._stop.is_set():
            rss = self._process.memory_info().rss / (1024 * 1024)
            if rss > self._peak:
                self._peak = rss
            self._stop.wait(self._interval)

    def stop(self):
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=2.0)
        rss = self._process.memory_info().rss / (1024 * 1024)
        if rss > self._peak:
            self._peak = rss

    @property
    def peak_mb(self):
        return self._peak

    @property
    def baseline_mb(self):
        return self._baseline

    @property
    def peak_delta_mb(self):
        return self._peak - self._baseline


def parse_expected_size(filename):
    inds_match = re.search(r'(\d+)inds', filename)
    snps_match = re.search(r'(\d+)([kKmM])?_snps', filename)
    inds = int(inds_match.group(1)) if inds_match else 0
    snps = 0
    if snps_match:
        val = int(snps_match.group(1))
        mult = snps_match.group(2)
        if mult:
            mult = mult.lower()
            if mult == 'k':
                snps = val * 1000
            elif mult == 'm':
                snps = val * 1000000
        else:
            snps = val
    return inds, snps


class BreakpointProvider:
    """Minimal object satisfying the `benchmark` interface expected by
    simulate_grg_recombination and simulate_numpy_recombination (needs
    .num_offspring_per_couple and .get_breakpoints())."""

    def __init__(self, num_offspring_per_couple):
        self.num_offspring_per_couple = num_offspring_per_couple

    def get_breakpoints(self, bp_range, expected_crossovers=1.5):
        num_bp = np.random.poisson(expected_crossovers)
        if num_bp == 0:
            return np.array([], dtype=int), num_bp
        low, high = bp_range
        length = high - low + 1
        num_bp = min(num_bp, length)
        k = num_bp * 3
        candidates = np.random.randint(low + 1, high, size=k)
        unique = np.unique(candidates)
        if unique.size < num_bp:
            extra_needed = num_bp - unique.size
            extra = np.random.randint(low + 1, high, size=extra_needed * 3)
            unique = np.unique(np.concatenate([unique, extra]))
        bp = np.sort(unique[:num_bp])
        return bp, num_bp


def system_info():
    import multiprocessing
    try:
        pygrgl_ver = pygrgl.__version__
    except AttributeError:
        pygrgl_ver = "unknown"
    return {
        "numpy_version": np.__version__,
        "pygrgl_version": pygrgl_ver,
        "platform": f"{platform.system()}-{platform.release()}-{platform.machine()}",
        "cpu_count": multiprocessing.cpu_count(),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "backend": _BACKEND,
    }


def run_grg(
    grg_file: Path,
    run_index: int,
    num_generations: int,
    offspring_per_couple: int,
    verification: bool = False,
    diagnostics: bool = False,
    serialize: bool = False,
    profile: bool = False,
):
    provider = BreakpointProvider(offspring_per_couple)
    fname = grg_file.name

    print(f"[GRG run {run_index}] Loading {fname}...")
    load_t0 = time.time()
    g = pygrgl.load_mutable_grg(str(grg_file))
    print(f"[GRG run {run_index}] Loaded in {time.time() - load_t0:.2f}s")

    base_samples = g.num_samples
    base_mutations = g.num_mutations
    base_nodes = g.num_nodes
    base_edges = g.num_edges
    base_genome = g.bp_range
    base_couples = base_samples // 2
    actual_offspring = base_couples * offspring_per_couple * num_generations

    result = {
        "file": fname,
        "method": "grg",
        "run_index": run_index,
        "num_samples_initial": base_samples,
        "num_snps": base_mutations,
        "base_nodes": base_nodes,
        "base_edges": base_edges,
        "bp_range": list(base_genome),
        "num_couples": base_couples,
        "num_offspring": actual_offspring,
        "num_generations": num_generations,
        "offspring_per_couple": offspring_per_couple,
    }

    structural_stats = None
    if diagnostics:
        print(f"[GRG run {run_index}] Computing structural stats...")
        t0 = time.perf_counter()
        structural_stats = compute_grg_structural_stats(g)
        print(f"[GRG run {run_index}] Structural stats done in {time.perf_counter() - t0:.2f}s")
        result["structural_stats"] = structural_stats
        del g
        gc.collect()
        g = pygrgl.load_mutable_grg(str(grg_file))

    # Measure pre-recombination graph size (off the timed path).
    # save_grg writes a compact serialization — this is the direct analog
    # to NumPy's matrix byte count.
    with tempfile.NamedTemporaryFile(suffix='.grg', delete=False) as _tf:
        _pre_grg_tmp = _tf.name
    try:
        pygrgl.save_grg(g, _pre_grg_tmp)
        grg_file_size_before_mb = os.path.getsize(_pre_grg_tmp) / (1024 * 1024)
    finally:
        if os.path.exists(_pre_grg_tmp):
            os.remove(_pre_grg_tmp)

    gc.collect()
    mem_before = get_process_memory_mb()
    peak_tracker = PeakRSSTracker(interval=0.1)
    peak_tracker.start()

    total_bp = 0
    all_offspring_ids = []
    liveness_snapshots = []
    offspring_ledger = [] if verification else None

    start = time.perf_counter()
    recomb = NonDuplicationRecombination(g)
    for gen in range(num_generations):
        print(f"  [Gen {gen+1}] Running GRG recombination...")
        offspring_ids, gen_bp = simulate_grg_recombination(
            provider, recomb, base_genome, N=base_genome[1],
            generation_index=gen, offspring_ledger=offspring_ledger,
        )
        all_offspring_ids.extend(offspring_ids)
        total_bp += gen_bp

        if diagnostics:
            t_pause = time.perf_counter()
            snap = compute_liveness(g)
            snap['gen'] = gen
            liveness_snapshots.append(snap)
            start += (time.perf_counter() - t_pause)

        if serialize:
            t_pause = time.perf_counter()
            pre_nodes = g.num_nodes
            pre_edges = g.num_edges
            with tempfile.NamedTemporaryFile(suffix='.grg', delete=False) as _tf:
                sr_tmp = _tf.name
            try:
                t0 = time.perf_counter()
                pygrgl.save_grg(g, sr_tmp)
                save_s = time.perf_counter() - t0
                file_size_mb = os.path.getsize(sr_tmp) / (1024 * 1024)
                t0 = time.perf_counter()
                g_reloaded = pygrgl.load_mutable_grg(sr_tmp)
                load_s = time.perf_counter() - t0
                post_nodes = g_reloaded.num_nodes
                post_edges = g_reloaded.num_edges
                result.setdefault("save_reload_per_gen", []).append({
                    "gen": gen,
                    "pre_nodes": pre_nodes,
                    "pre_edges": pre_edges,
                    "post_nodes": post_nodes,
                    "post_edges": post_edges,
                    "node_reduction_pct": 100.0 * (pre_nodes - post_nodes) / max(1, pre_nodes),
                    "edge_reduction_pct": 100.0 * (pre_edges - post_edges) / max(1, pre_edges),
                    "save_s": save_s,
                    "load_s": load_s,
                    "total_s": save_s + load_s,
                    "file_size_mb": file_size_mb,
                })
                del g_reloaded
            finally:
                if os.path.exists(sr_tmp):
                    os.remove(sr_tmp)
            gc.collect()
            start += time.perf_counter() - t_pause

        print(f"  [Gen {gen+1}] bp={gen_bp}, total_bp={total_bp}")

    elapsed = time.perf_counter() - start
    print(f"[GRG run {run_index}] Done in {elapsed:.4f}s")

    peak_tracker.stop()
    gc.collect()
    mem_after = get_process_memory_mb()

    # Measure post-recombination graph size (off the timed path).
    with tempfile.NamedTemporaryFile(suffix='.grg', delete=False) as _tf:
        _post_grg_tmp = _tf.name
    try:
        pygrgl.save_grg(g, _post_grg_tmp)
        grg_file_size_after_mb = os.path.getsize(_post_grg_tmp) / (1024 * 1024)
    finally:
        if os.path.exists(_post_grg_tmp):
            os.remove(_post_grg_tmp)

    result["time_ms"] = elapsed * 1000
    result["total_bp"] = total_bp
    result["mean_bp"] = total_bp / num_generations
    result["nodes_added"] = g.num_nodes - base_nodes
    result["edges_added"] = g.num_edges - base_edges
    result["memory_mb"] = mem_after
    result["memory_before_mb"] = mem_before
    result["memory_delta_mb"] = mem_after - mem_before
    result["memory_peak_mb"] = peak_tracker.peak_mb
    result["memory_peak_delta_mb"] = peak_tracker.peak_delta_mb
    result["grg_file_size_before_mb"] = grg_file_size_before_mb
    result["grg_file_size_after_mb"] = grg_file_size_after_mb
    result["grg_file_size_delta_mb"] = grg_file_size_after_mb - grg_file_size_before_mb

    print(f"[GRG run {run_index}] Memory: before={mem_before:.1f} MB, "
          f"after={mem_after:.1f} MB, delta={mem_after - mem_before:+.1f} MB, "
          f"peak={peak_tracker.peak_mb:.1f} MB, "
          f"peak_delta={peak_tracker.peak_delta_mb:+.1f} MB")
    print(f"[GRG run {run_index}] Graph size: before={grg_file_size_before_mb:.1f} MB, "
          f"after={grg_file_size_after_mb:.1f} MB, "
          f"delta={grg_file_size_after_mb - grg_file_size_before_mb:+.1f} MB")

    if diagnostics and liveness_snapshots:
        liveness_keys = (
            'total_nodes', 'alive', 'dead', 'dead_pct',
            'num_samples', 'dead_samples', 'dead_roots',
            'dead_with_mutations', 'dead_internal_empty', 'dead_mutation_count',
            'total_edges', 'alive_alive_edges', 'dead_alive_edges',
            'dead_dead_edges', 'dead_edges', 'dead_edges_pct',
        )
        n_snaps = len(liveness_snapshots)
        result["liveness"] = {
            "snapshot_count": n_snaps,
            "means": {
                k: sum(s[k] for s in liveness_snapshots) / n_snaps
                for k in liveness_keys
            },
            "snapshots": liveness_snapshots,
        }

    if verification:
        print(f"[GRG run {run_index}] Running verification...")
        audit_results = recomb.audit_check(raise_on_fail=False)
        audit_pass = sum(1 for r in audit_results.values() if r['pass'])
        print(f"  Audit: {audit_pass}/{len(audit_results)} pass")
        for r in audit_results.values():
            mark = "OK" if r['pass'] else "FAIL"
            print(f"    [{mark}] {r['desc']}: lhs={r['lhs']} rhs={r['rhs']}")

        mt_t0 = time.perf_counter()
        anc_count = compute_post_recomb_anc_counts(g)
        mt_result = check_offspring(g, all_offspring_ids, anc_count)
        mt_elapsed = time.perf_counter() - mt_t0
        mt_sum = mt_result['summary']
        print(f"  Multitree check: {mt_sum['violating']:,}/{mt_sum['total_checked']:,} "
              f"violate ({100*mt_sum['violation_rate']:.2f}%) in {mt_elapsed:.2f}s")

        ho_t0 = time.perf_counter()
        ho_result = check_offspring_haplotypes(g, offspring_ledger)
        ho_elapsed = time.perf_counter() - ho_t0
        print(f"  Haplotype oracle: {ho_result['num_pass']}/{ho_result['total_checked']} "
              f"pass in {ho_elapsed:.2f}s")
        if ho_result['num_fail'] > 0:
            print(f"  WARNING: {ho_result['num_fail']} offspring have incorrect mutations!")
            for fail in ho_result['failures'][:3]:
                print(f"    offspring {fail['offspring_id']}: "
                      f"{fail['missing_count']} missing, {fail['extra_count']} extra")

        result["verification"] = {
            "audit": audit_results,
            "multitree": {
                "wallclock_s": mt_elapsed,
                "offspring_checked": len(all_offspring_ids),
                "summary": mt_sum,
                "first_example": mt_result['first_example'],
            },
            "haplotype_oracle": {
                "wallclock_s": ho_elapsed,
                "total_checked": ho_result['total_checked'],
                "num_pass": ho_result['num_pass'],
                "num_fail": ho_result['num_fail'],
                "all_pass": ho_result['all_pass'],
                "failures": ho_result['failures'],
            },
        }

    if diagnostics:
        print(f"[GRG run {run_index}] Running instrumented diagnostic pass...")
        diag_g = pygrgl.load_mutable_grg(str(grg_file))
        diag_recomb = NonDuplicationRecombination(diag_g, instrument=True)
        init_caches_time_s = diag_recomb.stats["init_caches_time"]
        per_gen_stats = []
        diag_total_start = time.perf_counter()
        for gen in range(num_generations):
            diag_recomb.reset_stats()
            gen_start = time.perf_counter()
            simulate_grg_recombination(provider, diag_recomb, base_genome, N=base_genome[1])
            gen_elapsed = time.perf_counter() - gen_start
            snapshot = dict(diag_recomb.stats)
            snapshot["audit"] = dict(diag_recomb.audit)
            snapshot["generation_index"] = gen
            snapshot["generation_wallclock_s"] = gen_elapsed
            snapshot["num_nodes_after_gen"] = diag_g.num_nodes
            snapshot["num_edges_after_gen"] = diag_g.num_edges
            edges_added = (
                snapshot["audit"]["connect_calls_in_attach"]
                + snapshot["audit"]["connect_calls_in_extract"]
            )
            snapshot["edges_added_this_gen"] = edges_added
            per_gen_stats.append(snapshot)
            print(f"  [Diag gen {gen+1}] wallclock={gen_elapsed:.2f}s "
                  f"offspring={snapshot['offspring_count']} bubbles={snapshot['bubbles_created']}")
        diag_total = time.perf_counter() - diag_total_start

        aggregated_audit = NonDuplicationRecombination._fresh_audit()
        for snap in per_gen_stats:
            for k, v in snap.get("audit", {}).items():
                aggregated_audit[k] = aggregated_audit.get(k, 0) + v

        result["diagnostics"] = {
            "init_caches_time_s": init_caches_time_s,
            "diagnostic_total_wallclock_s": diag_total,
            "per_generation": per_gen_stats,
            "audit_aggregated": aggregated_audit,
        }
        del diag_g, diag_recomb
        gc.collect()

    if profile:
        import cProfile
        import pstats
        import io

        print(f"[GRG run {run_index}] Running cProfile pass...")
        prof_g = pygrgl.load_mutable_grg(str(grg_file))
        prof_recomb = NonDuplicationRecombination(prof_g)
        pr = cProfile.Profile()
        pr.enable()
        simulate_grg_recombination(provider, prof_recomb, base_genome, N=base_genome[1])
        pr.disable()

        stream = io.StringIO()
        ps = pstats.Stats(pr, stream=stream).strip_dirs()
        ps.sort_stats('cumtime').print_stats(30)
        cumtime_output = stream.getvalue()
        stream.seek(0); stream.truncate(0)
        ps.sort_stats('tottime').print_stats(30)
        tottime_output = stream.getvalue()

        print(cumtime_output)
        print(tottime_output)

        result["profile"] = {
            "cumtime_top30": cumtime_output,
            "tottime_top30": tottime_output,
        }
        del prof_g, prof_recomb
        gc.collect()

    del g, recomb
    gc.collect()
    return result


def run_numpy(
    grg_file: Path,
    run_index: int,
    num_generations: int,
    offspring_per_couple: int,
    memory_limit_mb: float = 15000.0,
):
    provider = BreakpointProvider(offspring_per_couple)
    fname = grg_file.name

    print(f"[NumPy run {run_index}] Loading {fname}...")
    g = pygrgl.load_mutable_grg(str(grg_file))

    base_samples = g.num_samples
    base_mutations = g.num_mutations
    base_genome = g.bp_range
    base_couples = base_samples // 2
    actual_offspring = base_couples * offspring_per_couple * num_generations

    est_mb = estimate_numpy_memory(base_mutations, base_samples + offspring_per_couple * base_couples)
    if est_mb > memory_limit_mb:
        print(f"[NumPy run {run_index}] Estimated memory ({est_mb:.1f} MB) exceeds limit "
              f"({memory_limit_mb} MB) — skipping.")
        return None

    result = {
        "file": fname,
        "method": "numpy",
        "run_index": run_index,
        "num_samples_initial": base_samples,
        "num_snps": base_mutations,
        "bp_range": list(base_genome),
        "num_couples": base_couples,
        "num_offspring": actual_offspring,
        "num_generations": num_generations,
        "offspring_per_couple": offspring_per_couple,
        "estimated_memory_mb": est_mb,
    }

    print(f"[NumPy run {run_index}] Building base population matrix...")
    base_matrix = grg_to_numpy_parallel(g)
    del g
    gc.collect()

    total_bp = 0
    offspring_matrix = base_matrix.copy()
    start = time.perf_counter()
    for gen in range(num_generations):
        print(f"  [Gen {gen+1}] Running NumPy recombination...")
        offspring_matrix, gen_bp = simulate_numpy_recombination(
            provider, offspring_matrix, bp_range=base_genome, expected_crossovers=1.5
        )
        total_bp += gen_bp
        print(f"  [Gen {gen+1}] bp={gen_bp}, total_bp={total_bp}")
    elapsed = time.perf_counter() - start
    print(f"[NumPy run {run_index}] Done in {elapsed:.4f}s")

    result["time_ms"] = elapsed * 1000
    result["total_bp"] = total_bp
    result["mean_bp"] = total_bp / num_generations
    result["memory_mb"] = est_mb
    result["nodes_added"] = 0
    result["edges_added"] = 0

    del base_matrix, offspring_matrix
    gc.collect()
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Execute a single benchmark run of one method on one .grg file."
    )
    parser.add_argument("--grg-file", type=Path, required=True,
                        help="Path to a single .grg file")
    parser.add_argument("--method", choices=["grg", "numpy"], required=True,
                        help="Which method to benchmark")
    parser.add_argument("--run-index", type=int, required=True,
                        help="Run index (0-based), used for output naming")
    parser.add_argument("--num-generations", type=int, default=1,
                        help="Sequential recombination generations (default: 1)")
    parser.add_argument("--offspring-per-couple", type=int, default=2,
                        help="Offspring per couple per generation (default: 2)")
    parser.add_argument("--output-dir", type=Path, default=Path("./output/runs"),
                        help="Directory for single-run result JSON (default: ./output/runs)")
    parser.add_argument("--memory-limit", type=float, default=15000.0,
                        help="Memory limit in MB for NumPy (default: 15000)")
    parser.add_argument("--verification", action="store_true",
                        help="Run audit + multitree verification (GRG only)")
    parser.add_argument("--diagnostics", action="store_true",
                        help="Run structural stats, liveness, and instrumented diag pass (GRG only)")
    parser.add_argument("--serialize", action="store_true",
                        help="Measure save_grg + load_mutable_grg per generation (GRG only)")
    parser.add_argument("--profile", action="store_true",
                        help="Run one generation under cProfile (GRG only)")
    parser.add_argument("--seed", type=int, default=None,
                        help="NumPy random seed for reproducibility")

    args = parser.parse_args()

    if args.seed is not None:
        np.random.seed(args.seed)

    if not args.grg_file.exists():
        print(f"Error: {args.grg_file} does not exist", file=sys.stderr)
        sys.exit(1)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    sysinfo = system_info()
    print(f"System: {sysinfo['platform']}, CPUs: {sysinfo['cpu_count']}, "
          f"Backend: {sysinfo['backend']}")

    if args.method == "grg":
        result = run_grg(
            grg_file=args.grg_file,
            run_index=args.run_index,
            num_generations=args.num_generations,
            offspring_per_couple=args.offspring_per_couple,
            verification=args.verification,
            diagnostics=args.diagnostics,
            serialize=args.serialize,
            profile=args.profile,
        )
    else:
        result = run_numpy(
            grg_file=args.grg_file,
            run_index=args.run_index,
            num_generations=args.num_generations,
            offspring_per_couple=args.offspring_per_couple,
            memory_limit_mb=args.memory_limit,
        )

    if result is None:
        print("Run skipped (memory limit exceeded).")
        sys.exit(0)

    result["system_info"] = sysinfo
    if args.seed is not None:
        result["seed"] = args.seed

    stem = args.grg_file.stem
    out_name = f"run_{stem}_{args.method}_{args.run_index}.json"
    out_path = args.output_dir / out_name
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2, default=str)

    print(f"Result written to {out_path}")


if __name__ == "__main__":
    main()
