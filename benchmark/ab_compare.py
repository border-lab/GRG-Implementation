#!/usr/bin/env python3
"""
ab_compare.py

Run NonDuplicationRecombination (interval-batched, current branch) and
NonDuplicationRecombinationLegacy (pre-change, copied from main) on the
same GRG file with identical parent shuffle + breakpoint draws, then
diff the instrumented stats so we can attribute differences purely to
the algorithmic change.

Usage:
    python3 ab_compare.py <grg_path> [--offspring-per-couple N]
                                     [--num-generations N]
                                     [--seed S]
                                     [--expected-crossovers X]
"""

import argparse
import gc
import json
import time
from pathlib import Path

import numpy as np
import pygrgl

from grg_recombination import (
    NonDuplicationRecombination,
    NonDuplicationRecombinationLegacy,
    recombination_intervals,
)


def get_breakpoints(bp_range, expected_crossovers):
    """Same logic as RecombinationBenchmarker.get_breakpoints; numpy-RNG-driven."""
    num_bp = np.random.poisson(expected_crossovers)
    if num_bp == 0:
        return np.array([], dtype=int), 0
    low, high = bp_range
    length = high - low + 1
    num_bp = min(num_bp, length)
    k = num_bp * 3
    candidates = np.random.randint(low + 1, high, size=k)
    unique = np.unique(candidates)
    if unique.size < num_bp:
        extra = np.random.randint(low + 1, high, size=(num_bp - unique.size) * 3)
        unique = np.unique(np.concatenate([unique, extra]))
    bp = np.sort(unique[:num_bp])
    return bp, num_bp


def run_one(klass, grg_path, seed, num_generations, num_offspring_per_couple,
            expected_crossovers):
    """Run one full simulation under `klass` with deterministic RNG."""
    np.random.seed(seed)
    g = pygrgl.load_mutable_grg(str(grg_path))
    base_nodes = g.num_nodes
    base_edges = g.num_edges
    base_genome = g.bp_range
    N = base_genome[1]

    recomb = klass(g, instrument=True)

    per_gen = []
    wall_total = 0.0
    for gen in range(num_generations):
        recomb.reset_stats()
        gen_start = time.perf_counter()

        # Replicate simulate_grg_recombination's deferral + shuffle pattern.
        samples = np.array(g.get_sample_nodes())
        np.random.shuffle(samples)
        new_offspring_ids = []
        prev_defer = recomb.defer_sample_updates
        recomb.defer_sample_updates = True
        try:
            for i in range(0, len(samples) - 1, 2):
                p1 = samples[i]
                p2 = samples[i + 1]
                for _ in range(num_offspring_per_couple):
                    bp, _num_bp = get_breakpoints(base_genome, expected_crossovers)
                    segments = recombination_intervals(p1, p2, bp, N)
                    off = recomb.recombine_multi(segments)
                    raw = recomb.NEGATIVE_NODE_IDS[abs(off) - 1]
                    new_offspring_ids.append(raw)
            new_offspring_ids.sort()
            g.set_samples(new_offspring_ids)
            recomb._pending_sample_removals.clear()
        finally:
            recomb.defer_sample_updates = prev_defer

        gen_elapsed = time.perf_counter() - gen_start
        wall_total += gen_elapsed
        snap = dict(recomb.stats)
        snap["generation_index"] = gen
        snap["generation_wallclock_s"] = gen_elapsed
        snap["num_nodes_after_gen"] = g.num_nodes
        snap["num_edges_after_gen"] = g.num_edges
        per_gen.append(snap)

    final = {
        "klass": klass.__name__,
        "wall_total_s": wall_total,
        "base_nodes": base_nodes,
        "base_edges": base_edges,
        "final_nodes": g.num_nodes,
        "final_edges": g.num_edges,
        "nodes_added": g.num_nodes - base_nodes,
        "edges_added": g.num_edges - base_edges,
        "per_generation": per_gen,
    }
    del recomb, g
    gc.collect()
    return final


def aggregate(per_gen):
    keys = ("recurse_attach_time", "apply_bubbles_time",
            "remove_mutation_time", "add_mutation_time", "get_mutation_by_id_time",
            "make_node_time", "connect_time",
            "recurse_attach_calls", "segments_processed",
            "bubbles_created", "mutations_moved", "visits_total",
            "make_node_calls", "connect_calls",
            "offspring_count")
    agg = {k: 0 for k in keys}
    for snap in per_gen:
        for k in keys:
            agg[k] += snap.get(k, 0)
    return agg


def fmt_pct(new, old):
    if old == 0:
        return "  n/a"
    return f"{(new - old) / old * 100:+6.2f}%"


def diff_report(legacy, new):
    print()
    print("=" * 90)
    print(f"A/B DIFF  ({legacy['klass']}  vs  {new['klass']})")
    print("=" * 90)

    L = aggregate(legacy["per_generation"])
    N = aggregate(new["per_generation"])

    rows = [
        ("wallclock_total_s",   legacy["wall_total_s"],          new["wall_total_s"]),
        ("nodes_added",         legacy["nodes_added"],           new["nodes_added"]),
        ("edges_added",         legacy["edges_added"],           new["edges_added"]),
        ("",                    "",                              ""),
        ("recurse_attach_calls", L["recurse_attach_calls"],      N["recurse_attach_calls"]),
        ("segments_processed",  L["segments_processed"],         N["segments_processed"]),
        ("bubbles_created",     L["bubbles_created"],            N["bubbles_created"]),
        ("mutations_moved",     L["mutations_moved"],            N["mutations_moved"]),
        ("visits_total",        L["visits_total"],               N["visits_total"]),
        ("make_node_calls",     L["make_node_calls"],            N["make_node_calls"]),
        ("connect_calls",       L["connect_calls"],              N["connect_calls"]),
        ("",                    "",                              ""),
        ("recurse_attach_time", L["recurse_attach_time"],        N["recurse_attach_time"]),
        ("apply_bubbles_time",  L["apply_bubbles_time"],         N["apply_bubbles_time"]),
        ("remove_mutation_time", L["remove_mutation_time"],      N["remove_mutation_time"]),
        ("add_mutation_time",   L["add_mutation_time"],          N["add_mutation_time"]),
        ("get_mutation_by_id_time", L["get_mutation_by_id_time"], N["get_mutation_by_id_time"]),
    ]

    print(f"{'metric':28s}  {'legacy':>14s}  {'new':>14s}  {'delta':>12s}")
    print("-" * 90)
    for name, lv, nv in rows:
        if name == "":
            print()
            continue
        if isinstance(lv, float) or isinstance(nv, float):
            print(f"{name:28s}  {lv:14.4f}  {nv:14.4f}  {fmt_pct(nv, lv):>12s}")
        else:
            print(f"{name:28s}  {lv:14d}  {nv:14d}  {fmt_pct(nv, lv):>12s}")

    # Per-offspring sanity checks (offspring count must match for fair A/B).
    if L["offspring_count"] != N["offspring_count"]:
        print(f"\nWARNING: offspring counts differ! legacy={L['offspring_count']} new={N['offspring_count']}")
    else:
        print(f"\noffspring per implementation: {L['offspring_count']}")
        print(f"bubbles/offspring  legacy={L['bubbles_created']/max(1,L['offspring_count']):.3f}  "
              f"new={N['bubbles_created']/max(1,N['offspring_count']):.3f}")
        print(f"muts_moved/offspring  legacy={L['mutations_moved']/max(1,L['offspring_count']):.2f}  "
              f"new={N['mutations_moved']/max(1,N['offspring_count']):.2f}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("grg_path", type=Path)
    p.add_argument("--offspring-per-couple", type=int, default=2)
    p.add_argument("--num-generations", type=int, default=1)
    p.add_argument("--seed", type=int, default=20260502)
    p.add_argument("--expected-crossovers", type=float, default=1.5)
    p.add_argument("--out", type=Path, default=None,
                   help="Optional path to write JSON dump of both runs.")
    args = p.parse_args()

    if not args.grg_path.exists():
        raise SystemExit(f"GRG file not found: {args.grg_path}")

    # File-level info.
    g = pygrgl.load_mutable_grg(str(args.grg_path))
    print(f"file: {args.grg_path.name}")
    print(f"  nodes={g.num_nodes}  edges={g.num_edges}  samples={g.num_samples}  "
          f"mutations={g.num_mutations}  bp_range={g.bp_range}")
    del g
    gc.collect()

    common = dict(
        grg_path=args.grg_path,
        seed=args.seed,
        num_generations=args.num_generations,
        num_offspring_per_couple=args.offspring_per_couple,
        expected_crossovers=args.expected_crossovers,
    )

    print(f"\n[legacy] running...")
    legacy = run_one(NonDuplicationRecombinationLegacy, **common)
    print(f"  wall={legacy['wall_total_s']:.3f}s  nodes_added={legacy['nodes_added']}")

    print(f"\n[new]    running...")
    new = run_one(NonDuplicationRecombination, **common)
    print(f"  wall={new['wall_total_s']:.3f}s  nodes_added={new['nodes_added']}")

    diff_report(legacy, new)

    if args.out:
        with open(args.out, "w") as f:
            json.dump({"legacy": legacy, "new": new, "args": vars(args) | {"grg_path": str(args.grg_path)}}, f, indent=2, default=str)
        print(f"\nwrote: {args.out}")


if __name__ == "__main__":
    main()
