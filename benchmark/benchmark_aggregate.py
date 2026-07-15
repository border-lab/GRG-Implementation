#!/usr/bin/env python3
"""
Benchmark result aggregator.

Reads single-run JSON files produced by benchmark_run.py, groups them by
(file, method), computes summary statistics (mean/std/min/max), and writes
a combined CSV + JSON that matches the format of the original
benchmark_recombination.py output.

Example
-------
python benchmark_aggregate.py \\
    --input-dir ./output/runs \\
    --output-dir ./output \\
    --output-prefix benchmark_recombination_results
"""

import os
import json
import csv
import argparse
import numpy as np
from pathlib import Path
from collections import defaultdict


CSV_COLUMNS = [
    "file",
    "num_offspring",
    "implementation",
    "num_samples_initial",
    "num_snps",
    "num_runs",
    "num_bp",
    "mean_bp",
    "mean_time_ms",
    "std_time_ms",
    "min_time_ms",
    "max_time_ms",
    "nodes_added",
    "edges_added",
    "memory_mb",
    "mean_dead_nodes",
    "mean_dead_pct",
    "mean_dead_mutations",
    "mean_dead_edges",
    "mean_dead_edges_pct",
]


def load_run_files(input_dir: Path):
    files = sorted(input_dir.glob("run_*.json"))
    if not files:
        print(f"No run_*.json files found in {input_dir}")
        return []
    runs = []
    for f in files:
        with open(f) as fh:
            runs.append(json.load(fh))
    print(f"Loaded {len(runs)} run files from {input_dir}")
    return runs


def group_runs(runs):
    groups = defaultdict(list)
    for r in runs:
        key = (r["file"], r["method"])
        groups[key].append(r)
    return groups


def aggregate_grg(runs, implementation="GRG Native"):
    times = [r["time_ms"] for r in runs]
    total_bp = sum(r["total_bp"] for r in runs)
    num_runs = len(runs)
    first = runs[0]
    num_generations = first["num_generations"]

    nodes_added_avg = np.mean([r["nodes_added"] for r in runs])
    edges_added_avg = np.mean([r["edges_added"] for r in runs])
    memory_avg = np.mean([r["memory_mb"] for r in runs])

    liveness_means = {
        "dead": 0.0,
        "dead_pct": 0.0,
        "dead_mutation_count": 0.0,
        "dead_edges": 0.0,
        "dead_edges_pct": 0.0,
    }
    all_liveness = []
    for r in runs:
        liv = r.get("liveness")
        if liv and liv.get("snapshots"):
            all_liveness.extend(liv["snapshots"])
    if all_liveness:
        n = len(all_liveness)
        for k in liveness_means:
            liveness_means[k] = sum(s.get(k, 0) for s in all_liveness) / n

    row = {
        "file": first["file"],
        "num_offspring": first["num_offspring"],
        "implementation": implementation,
        "num_samples_initial": first["num_samples_initial"],
        "num_snps": first["num_snps"],
        "num_runs": num_runs,
        "num_bp": total_bp,
        "mean_bp": total_bp / (num_runs * num_generations),
        "mean_time_ms": float(np.mean(times)),
        "std_time_ms": float(np.std(times)),
        "min_time_ms": float(np.min(times)),
        "max_time_ms": float(np.max(times)),
        "nodes_added": nodes_added_avg,
        "edges_added": edges_added_avg,
        "memory_mb": memory_avg,
        "mean_dead_nodes": liveness_means["dead"],
        "mean_dead_pct": liveness_means["dead_pct"],
        "mean_dead_mutations": liveness_means["dead_mutation_count"],
        "mean_dead_edges": liveness_means["dead_edges"],
        "mean_dead_edges_pct": liveness_means["dead_edges_pct"],
    }
    return row


def aggregate_numpy(runs):
    times = [r["time_ms"] for r in runs]
    total_bp = sum(r["total_bp"] for r in runs)
    num_runs = len(runs)
    first = runs[0]
    num_generations = first["num_generations"]

    row = {
        "file": first["file"],
        "num_offspring": first["num_offspring"],
        "implementation": "NumPy Baseline",
        "num_samples_initial": first["num_samples_initial"],
        "num_snps": first["num_snps"],
        "num_runs": num_runs,
        "num_bp": total_bp,
        "mean_bp": total_bp / (num_runs * num_generations),
        "mean_time_ms": float(np.mean(times)),
        "std_time_ms": float(np.std(times)),
        "min_time_ms": float(np.min(times)),
        "max_time_ms": float(np.max(times)),
        "nodes_added": 0,
        "edges_added": 0,
        "memory_mb": first.get("estimated_memory_mb", first.get("memory_mb", 0.0)),
        "mean_dead_nodes": 0.0,
        "mean_dead_pct": 0.0,
        "mean_dead_mutations": 0.0,
        "mean_dead_edges": 0.0,
        "mean_dead_edges_pct": 0.0,
    }
    return row


def collect_diagnostics(groups):
    diagnostics = {}
    for (fname, method), runs in groups.items():
        if method not in ("grg", "grg_parallel"):
            continue
        file_diag = {}
        for r in runs:
            if "diagnostics" in r:
                file_diag["per_generation"] = r["diagnostics"].get("per_generation")
                file_diag["audit_aggregated"] = r["diagnostics"].get("audit_aggregated")
                file_diag["diagnostic_total_wallclock_s"] = r["diagnostics"].get("diagnostic_total_wallclock_s")
                file_diag["init_caches_time_s"] = r["diagnostics"].get("init_caches_time_s")
            if "structural_stats" in r:
                file_diag["structural"] = r["structural_stats"]
            if "verification" in r:
                file_diag.setdefault("verification_checks", []).append({
                    "run_index": r["run_index"],
                    "audit_identities": r["verification"]["audit"],
                    "multitree": r["verification"]["multitree"],
                })
            if "save_reload_per_gen" in r:
                file_diag["save_reload_per_gen"] = r["save_reload_per_gen"]
            if "liveness" in r:
                file_diag.setdefault("liveness_snapshots_all_runs", []).extend(
                    r["liveness"].get("snapshots", [])
                )
            if "profile" in r:
                file_diag["profile"] = r["profile"]
        if file_diag:
            # Keep the plain filename as key for "grg" (backward compatible
            # with existing JSON consumers); disambiguate grg_parallel so it
            # doesn't clobber the sequential entry when both ran on the same file.
            key = fname if method == "grg" else f"{fname} [{method}]"
            diagnostics[key] = file_diag
    return diagnostics


def build_config(runs):
    first = runs[0]
    return {
        "timed_runs": len([r for r in runs if r["file"] == first["file"] and r["method"] == first["method"]]),
        "num_generations": first["num_generations"],
        "num_offspring_per_couple": first["offspring_per_couple"],
        "num_offspring": first["num_offspring"],
    }


def save_results(rows, groups, diagnostics, runs, output_dir, prefix):
    output_dir.mkdir(parents=True, exist_ok=True)

    csv_path = output_dir / f"{prefix}.csv"
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    sysinfo = {}
    for r in runs:
        if "system_info" in r:
            sysinfo = r["system_info"]
            break

    json_path = output_dir / f"{prefix}.json"
    with open(json_path, 'w') as f:
        json.dump({
            "system_info": sysinfo,
            "config": build_config(runs),
            "results": rows,
            "diagnostics": diagnostics,
        }, f, indent=2, default=str)

    print(f"Results saved to {csv_path} and {json_path}")


def print_summary(rows):
    print("\n" + "=" * 80)
    print("Aggregated Results")
    print("=" * 80)
    for row in rows:
        impl = row["implementation"]
        fname = row["file"]
        n = row["num_runs"]
        mean = row["mean_time_ms"]
        std = row["std_time_ms"]
        print(f"  {fname} | {impl:15s} | {n} runs | {mean:.2f}ms +/- {std:.2f}ms")

    grg_rows = {r["file"]: r for r in rows if r["implementation"] == "GRG Native"}
    np_rows = {r["file"]: r for r in rows if r["implementation"] == "NumPy Baseline"}
    for fname in grg_rows:
        if fname in np_rows:
            speedup = np_rows[fname]["mean_time_ms"] / grg_rows[fname]["mean_time_ms"]
            print(f"  {fname} | Speedup: {speedup:.2f}x (NumPy / GRG)")
    print("=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate single-run benchmark JSONs into CSV + JSON summary."
    )
    parser.add_argument("--input-dir", type=Path, required=True,
                        help="Directory containing run_*.json files from benchmark_run.py")
    parser.add_argument("--output-dir", type=Path, default=Path("."),
                        help="Directory to write aggregated CSV and JSON (default: .)")
    parser.add_argument("--output-prefix", type=str,
                        default="benchmark_recombination_results",
                        help="Filename prefix for output files (default: benchmark_recombination_results)")

    args = parser.parse_args()

    runs = load_run_files(args.input_dir)
    if not runs:
        return

    groups = group_runs(runs)

    file_order = []
    seen = set()
    for r in runs:
        if r["file"] not in seen:
            file_order.append(r["file"])
            seen.add(r["file"])

    rows = []
    for fname in file_order:
        grg_key = (fname, "grg")
        grg_parallel_key = (fname, "grg_parallel")
        np_key = (fname, "numpy")
        if grg_key in groups:
            rows.append(aggregate_grg(groups[grg_key], implementation="GRG Native"))
        if grg_parallel_key in groups:
            rows.append(aggregate_grg(groups[grg_parallel_key], implementation="GRG Parallel"))
        if np_key in groups:
            rows.append(aggregate_numpy(groups[np_key]))

    diagnostics = collect_diagnostics(groups)

    print_summary(rows)
    save_results(rows, groups, diagnostics, runs, args.output_dir, args.output_prefix)


if __name__ == "__main__":
    main()
