"""
Diagnostic profiler for the GRG recombination pathology investigation.

Compares the structural properties of two GRG files (e.g. 500inds_1M vs
5000inds_1M) and measures per-offspring `|V_visited|` plus a cProfile
breakdown of the recombination hot path.

Usage:
    python profile_pathology.py [files...]   (default: 500inds_1M and 5000inds_1M)
"""

import sys
import time
import cProfile
import pstats
import io
from pathlib import Path
from collections import deque

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from pygrgl import load_mutable_grg
from grg_recombination import NonDuplicationRecombination, recombination_intervals


def graph_stats(grg):
    """Quick structural metrics. Avoids touching every node twice."""
    n = grg.num_nodes
    samples = [i for i in range(n) if grg.is_sample(i)]
    roots = list(grg.get_root_nodes())

    # Topological depth: BFS down from roots, tracking longest path.
    depth = {r: 0 for r in roots}
    queue = deque(roots)
    while queue:
        node = queue.popleft()
        d = depth[node]
        for child in grg.get_down_edges(node):
            cd = depth.get(child, -1)
            if cd < d + 1:
                depth[child] = d + 1
                queue.append(child)

    sample_depths = [depth.get(s, 0) for s in samples]

    # Fanouts: sample to keep this cheap.
    up_total = 0
    down_total = 0
    up_max = 0
    for i in range(n):
        u = len(grg.get_up_edges(i))
        d = len(grg.get_down_edges(i))
        up_total += u
        down_total += d
        if u > up_max:
            up_max = u

    return {
        "num_nodes": n,
        "num_edges": grg.num_edges,
        "num_samples": len(samples),
        "num_mutations": grg.num_mutations,
        "num_roots": len(roots),
        "internal_per_sample": (n - len(samples)) / max(1, len(samples)),
        "avg_up_fanout": up_total / max(1, n),
        "max_up_fanout": up_max,
        "avg_down_fanout": down_total / max(1, n),
        "sample_depth_min": min(sample_depths) if sample_depths else 0,
        "sample_depth_avg": sum(sample_depths) / max(1, len(sample_depths)),
        "sample_depth_max": max(sample_depths) if sample_depths else 0,
    }


def random_breakpoints(N, expected_crossovers=1.5):
    k = max(1, int(np.random.poisson(expected_crossovers)))
    bps = sorted(set(int(x) for x in np.random.randint(1, N, k)))
    return bps


def profile_file(path, num_offspring=20, expected_crossovers=1.5):
    print(f"\n{'=' * 70}")
    print(f"  {path.name}")
    print(f"{'=' * 70}")

    t = time.perf_counter()
    grg = load_mutable_grg(str(path))
    print(f"loaded in {time.perf_counter() - t:.2f}s")

    t = time.perf_counter()
    stats = graph_stats(grg)
    print(f"\nstructural stats (computed in {time.perf_counter() - t:.2f}s):")
    for k, v in stats.items():
        if isinstance(v, float):
            print(f"  {k:30s} {v:,.2f}")
        else:
            print(f"  {k:30s} {v:,}")

    # Set up recombiner and instrument recombine_multi to count visits.
    recomb = NonDuplicationRecombination(grg)
    recomb.defer_sample_updates = True

    samples = list(grg.get_sample_nodes())
    np.random.shuffle(samples)
    if len(samples) < 2:
        print("not enough samples")
        return

    visit_counts = []
    bubble_counts = []
    timings = []

    original_recombine = recomb.recombine_multi

    def counted_recombine(segments):
        gen_before = recomb._gen_visited
        bubbles_before = len(recomb._pending_bubbles)
        t0 = time.perf_counter()
        result = original_recombine(segments)
        timings.append(time.perf_counter() - t0)
        gen_after = recomb._gen_visited
        # _recurse_attach bumps _gen_visited once per call. Count nodes
        # whose _visited_gen entry falls in (gen_before, gen_after].
        gen_set = set(range(gen_before + 1, gen_after + 1))
        count = 0
        for v in recomb._visited_gen:
            if v in gen_set:
                count += 1
        visit_counts.append(count)
        # bubbles_added in recombine_multi: pending_bubbles is cleared
        # at start, so size after-traversal-before-apply is the count.
        # We need to wrap differently to capture it; cheat by inspecting
        # _modified_nodes delta won't work either (cleared at end).
        # Approximate via NEGATIVE_NODE_IDS growth + bubble count from
        # this run -- skip for simplicity.
        return result

    recomb.recombine_multi = counted_recombine

    bp_range = grg.bp_range
    N = bp_range[1]
    nodes_before = grg.num_nodes

    print(f"\nrunning {num_offspring} offspring with cProfile...")
    profiler = cProfile.Profile()
    profiler.enable()

    for i in range(num_offspring):
        p1 = int(samples[(2 * i) % len(samples)])
        p2 = int(samples[(2 * i + 1) % len(samples)])
        bps = random_breakpoints(N, expected_crossovers)
        segments = recombination_intervals(p1, p2, bps, N)
        recomb.recombine_multi(segments)

    profiler.disable()
    nodes_after = grg.num_nodes

    visits = visit_counts
    print(f"\n|V_visited| per offspring:")
    print(f"  min:   {min(visits):,}")
    print(f"  avg:   {sum(visits) / len(visits):,.0f}")
    print(f"  max:   {max(visits):,}")
    print(f"  total: {sum(visits):,}")
    print(f"\nper-offspring wallclock:")
    print(f"  min:   {min(timings) * 1000:.2f} ms")
    print(f"  avg:   {sum(timings) * 1000 / len(timings):.2f} ms")
    print(f"  max:   {max(timings) * 1000:.2f} ms")
    print(f"\nnodes added: {nodes_after - nodes_before:,} ({(nodes_after - nodes_before) / num_offspring:.1f} per offspring)")

    # Per-visit cost (microseconds)
    total_seconds = sum(timings)
    total_visits = sum(visits)
    if total_visits:
        print(f"per-visit cost: {total_seconds * 1e6 / total_visits:.2f} us")

    print(f"\ntop functions by cumulative time:")
    s = io.StringIO()
    ps = pstats.Stats(profiler, stream=s).sort_stats("cumulative")
    ps.print_stats(15)
    print(s.getvalue())


def main():
    base = Path(__file__).parent / "test_grg_files"
    if len(sys.argv) > 1:
        files = [Path(p) for p in sys.argv[1:]]
    else:
        files = [
            base / "500inds_1M_snps.grg",
            base / "5000inds_1M_snps.grg",
        ]

    for f in files:
        if not f.exists():
            print(f"missing: {f}")
            continue
        profile_file(f, num_offspring=20)


if __name__ == "__main__":
    main()
