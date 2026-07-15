"""
multitree_check.py

Approach B from the multitree-verification discussion: per-offspring cardinality
test, focused on whether one generation of recombination introduces diamond
patterns. Assumes the input GRG is already a multitree.

For each offspring o, the test is:

    expected[o] = sum over up-edges p of (1 + anc_count[p])
    actual[o]   = | unique ancestors of o |

In a multitree, expected == actual. If expected > actual, two of o's parents
share an ancestor and the difference is the count of redundant paths.

`anc_count[v]` is computed by DP on the *post-recomb* GRG iterating high-to-low
node ID. Correctness rests on three observations:
  - Bubbles have no up-edges, so anc_count[bubble] = 0 (processed first,
    highest IDs).
  - Pre-existing nodes that gained a bubble as a new up-edge during recomb
    don't violate multitree (bubbles introduce no overlap), so the DP gives
    the correct post-recomb true ancestor count for them.
  - Offspring nodes' DP values come out wrong because their pre-existing-node
    parents have lower IDs and aren't yet processed. But offspring are leaves
    -- nothing references their anc_count -- so the error doesn't propagate.
    The cardinality test for each offspring uses only its parents'
    (already-correct) anc_count values.

Run as:
    python multitree_check.py path/to/file.grg [--seed 42] \\
            [--offspring-per-couple 2] [--output-json out.json]
"""

import argparse
import json
import sys
import time
from collections import deque
from pathlib import Path

import numpy as np
import pygrgl


class BenchmarkStub:
    """Minimal stand-in for the benchmark instance that
    `simulate_grg_recombination` expects. Mirrors
    `RecombinationBenchmarker.get_breakpoints` from
    benchmark_recombination.py."""

    def __init__(self, num_offspring_per_couple: int = 2):
        self.num_offspring_per_couple = num_offspring_per_couple

    def get_breakpoints(self, bp_range, expected_crossovers: float = 1.5):
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


def compute_post_recomb_anc_counts(grg) -> list:
    """Per-node true ancestor count on the post-recomb GRG via Kahn's
    topological-order DP.

    For each node v in topological order (parents-before-children):
        anc_count[v] = Σ_p (1 + anc_count[p])  over up-edges p

    This equals |ancestors(v)| iff the *post-recomb* DAG is a multitree
    (parents' ancestor sets are disjoint); otherwise it overcounts, and the
    overcount is exactly the number of redundant paths in the diamond.

    Why topo-order, not ID-order: pre-recomb GRGs have parents at higher IDs,
    so iterating high-to-low works. Post-recomb, direct-attach and
    path-compression-attach add edges where the offspring (high ID) has a
    pre-existing parent (low ID). In multi-generation runs, the gen-1
    offspring then becomes a parent of a gen-2 offspring, so gen-1 offspring
    anc_counts get read -- and the ID-based DP undercounts them whenever
    they have low-ID pre-existing parents.

    Implementation: Kahn's algorithm seeded from true roots (no up-edges),
    propagating along down-edges. O(|V| + |E|)."""
    n = grg.num_nodes
    anc_count = [0] * n
    indeg = [0] * n
    for v in range(n):
        indeg[v] = len(grg.get_up_edges(v))

    queue = deque(v for v in range(n) if indeg[v] == 0)

    while queue:
        v = queue.popleft()
        total = 0
        for p in grg.get_up_edges(v):
            total += 1 + anc_count[p]
        anc_count[v] = total
        for c in grg.get_down_edges(v):
            indeg[c] -= 1
            if indeg[c] == 0:
                queue.append(c)
    return anc_count


def bfs_unique_ancestor_count(grg, start: int) -> int:
    """Number of unique ancestors of `start`, walking up via `get_up_edges`."""
    visited = {start}
    queue = deque([start])
    while queue:
        n = queue.popleft()
        for p in grg.get_up_edges(n):
            if p not in visited:
                visited.add(p)
                queue.append(p)
    return len(visited) - 1


def bfs_collect_ancestors(grg, start: int) -> set:
    """Same as bfs_unique_ancestor_count but returns the ancestor set itself.
    Used to characterize a single violation example; not on the hot path."""
    visited = {start}
    queue = deque([start])
    while queue:
        n = queue.popleft()
        for p in grg.get_up_edges(n):
            if p not in visited:
                visited.add(p)
                queue.append(p)
    visited.discard(start)
    return visited


def compute_liveness(grg, samples=None, return_visited: bool = False) -> dict:
    """Count nodes reachable upward from `samples` ("alive") vs not ("dead").

    A node is alive if and only if it is an ancestor of (or equal to) at least
    one node in `samples`. Anything else is "deadweight" the GRG carries
    forward without contributing to any current sample's haplotype.

    Single shared bytearray visited mask, all samples seeded at once -- one
    BFS over the up-edges, not one per sample. O(|V| + |E|).

    Args:
        grg: pygrgl.MutableGRG (or any object exposing num_nodes,
            get_up_edges, get_sample_nodes, is_sample, get_mutations_for_node).
        samples: iterable of node IDs to seed from. Defaults to the GRG's
            current sample set. Pass an explicit set to measure liveness
            against arbitrary seeds (e.g. founder lineages) without touching
            set_samples.
        return_visited: if True, also include the visited bytearray in the
            result under key 'visited' (useful for assertions in tests).

    Returns a dict with: total_nodes, alive, dead, dead_pct, num_samples,
    dead_samples, dead_roots, dead_with_mutations, dead_internal_empty,
    dead_mutation_count, total_edges, alive_alive_edges, dead_alive_edges,
    dead_dead_edges, dead_edges, dead_edges_pct.

    Edge classification: every (parent, child) edge in the GRG falls into
    one of three buckets based on endpoint liveness:
      - alive_alive_edges: both endpoints reachable from samples (useful).
      - dead_alive_edges:  child dead, parent alive (dangling -- the live
                           parent still has a down-edge to a stranded child).
      - dead_dead_edges:   both endpoints dead (entirely in a dead subgraph).
    The fourth combination (child alive, parent dead) is impossible: BFS-
    upward from samples reaches every ancestor of every alive node, so an
    alive child always has all-alive parents.

    dead_samples must always be 0 (a sample reaches itself trivially); a
    non-zero value indicates a malformed input."""
    if samples is None:
        samples = grg.get_sample_nodes()
    samples = list(samples)

    n = grg.num_nodes
    visited = bytearray(n)
    queue = deque()
    for s in samples:
        if 0 <= s < n and not visited[s]:
            visited[s] = 1
            queue.append(s)

    while queue:
        node = queue.popleft()
        for p in grg.get_up_edges(node):
            if not visited[p]:
                visited[p] = 1
                queue.append(p)

    alive = sum(visited)
    dead = n - alive

    dead_samples = 0
    dead_roots = 0
    dead_with_mutations = 0
    dead_internal_empty = 0
    dead_mutation_count = 0
    alive_alive_edges = 0
    dead_alive_edges = 0
    dead_dead_edges = 0

    for node in range(n):
        parents = grg.get_up_edges(node)
        if visited[node]:
            # By the BFS-upward invariant, every parent of an alive node is
            # also alive -- no per-parent visit needed.
            alive_alive_edges += len(parents)
            continue
        # Dead node: classify content (sample/root/mut-bearing/empty) and
        # bucket each up-edge by parent liveness.
        if grg.is_sample(node):
            dead_samples += 1
        if not parents:
            dead_roots += 1
        muts = grg.get_mutations_for_node(node)
        num_muts = len(muts)
        if num_muts:
            dead_with_mutations += 1
            dead_mutation_count += num_muts
        else:
            dead_internal_empty += 1
        for p in parents:
            if visited[p]:
                dead_alive_edges += 1
            else:
                dead_dead_edges += 1

    total_edges = alive_alive_edges + dead_alive_edges + dead_dead_edges
    dead_edges = dead_alive_edges + dead_dead_edges

    result = {
        'total_nodes':         n,
        'alive':               alive,
        'dead':                dead,
        'dead_pct':            (100.0 * dead / n) if n else 0.0,
        'num_samples':         len(samples),
        'dead_samples':        dead_samples,
        'dead_roots':          dead_roots,
        'dead_with_mutations': dead_with_mutations,
        'dead_internal_empty': dead_internal_empty,
        'dead_mutation_count': dead_mutation_count,
        'total_edges':         total_edges,
        'alive_alive_edges':   alive_alive_edges,
        'dead_alive_edges':    dead_alive_edges,
        'dead_dead_edges':     dead_dead_edges,
        'dead_edges':          dead_edges,
        'dead_edges_pct':      (100.0 * dead_edges / total_edges) if total_edges else 0.0,
    }
    if return_visited:
        result['visited'] = visited
    return result


def find_overlapping_ancestors(grg, parents: list) -> dict:
    """For one violation example: collect each parent's ancestor set (including
    the parent itself) and find ancestors reached via more than one parent."""
    parent_sets = []
    for p in parents:
        anc = bfs_collect_ancestors(grg, p)
        anc.add(p)
        parent_sets.append(anc)

    counts = {}
    for s in parent_sets:
        for a in s:
            counts[a] = counts.get(a, 0) + 1
    overlapping = {a: c for a, c in counts.items() if c > 1}
    return {
        'overlapping': overlapping,
        'parent_anc_sizes': [len(s) for s in parent_sets],
    }


def check_offspring(grg, offspring_ids, anc_count) -> dict:
    """Run the cardinality test for every offspring; return summary + first
    violation example (if any).

    `anc_count` is the post-recomb DP table from `compute_post_recomb_anc_counts`.
    Indexes for pre-existing parents and bubble parents are both correct;
    offspring's own anc_count slot is wrong but never read."""
    violations = []
    skipped_no_parents = 0

    for o in offspring_ids:
        parents = list(grg.get_up_edges(o))
        if not parents:
            skipped_no_parents += 1
            continue

        expected = sum(1 + anc_count[p] for p in parents)
        actual = bfs_unique_ancestor_count(grg, o)

        if expected != actual:
            violations.append({
                'offspring': int(o),
                'num_parents': len(parents),
                'parents': [int(p) for p in parents],
                'expected_ancestors': int(expected),
                'actual_ancestors': int(actual),
                'extra_paths': int(expected - actual),
            })

    extras = [v['extra_paths'] for v in violations]
    summary = {
        'total_checked': len(offspring_ids) - skipped_no_parents,
        'skipped_no_parents': skipped_no_parents,
        'violating': len(violations),
        'violation_rate': (len(violations) / max(1, len(offspring_ids) - skipped_no_parents)),
        'total_redundant_paths': int(sum(extras)),
    }
    if extras:
        sorted_extras = sorted(extras)
        m = len(sorted_extras)
        summary['extra_paths'] = {
            'min': sorted_extras[0],
            'p50': sorted_extras[m // 2],
            'p95': sorted_extras[min(m - 1, int(m * 0.95))],
            'p99': sorted_extras[min(m - 1, int(m * 0.99))],
            'max': sorted_extras[-1],
            'mean': sum(sorted_extras) / m,
        }

    first_example = None
    if violations:
        v = violations[0]
        overlap_info = find_overlapping_ancestors(grg, v['parents'])
        # Cap reported overlapping IDs so output stays readable on big diamonds.
        sample_overlap = sorted(overlap_info['overlapping'].keys())[:20]
        first_example = {
            **v,
            'parent_anc_sizes': overlap_info['parent_anc_sizes'],
            'num_overlapping_ancestors': len(overlap_info['overlapping']),
            'sample_overlapping_ids': [int(a) for a in sample_overlap],
        }

    return {
        'summary': summary,
        'first_example': first_example,
        'violations': violations,
    }


def print_report(grg_path: Path, initial_num_nodes: int, num_offspring: int,
                 result: dict, runtime_s: float):
    s = result['summary']
    print(f"\nMultitree check on {grg_path.name}:")
    print(f"  Initial nodes: {initial_num_nodes:,}")
    print(f"  Offspring checked: {s['total_checked']:,}")
    if s['skipped_no_parents']:
        print(f"  Skipped (no parents): {s['skipped_no_parents']:,}")
    print(f"  Violations: {s['violating']:,} offspring "
          f"({100*s['violation_rate']:.2f}%)")
    if 'extra_paths' in s:
        ep = s['extra_paths']
        print(f"  Extra-paths per violating offspring:")
        print(f"    min={ep['min']}, p50={ep['p50']}, p95={ep['p95']}, "
              f"p99={ep['p99']}, max={ep['max']}, mean={ep['mean']:.2f}")
    print(f"  Total redundant paths across all offspring: "
          f"{s['total_redundant_paths']:,}")

    fe = result['first_example']
    if fe is not None:
        print(f"\n  Sample violation: offspring {fe['offspring']}")
        print(f"    {fe['num_parents']} parents")
        print(f"    expected ancestors: {fe['expected_ancestors']:,}")
        print(f"    actual ancestors:   {fe['actual_ancestors']:,}")
        print(f"    extra paths: {fe['extra_paths']:,}")
        print(f"    overlapping ancestors (count): "
              f"{fe['num_overlapping_ancestors']:,}")
        if fe['sample_overlapping_ids']:
            print(f"    overlapping ancestor IDs (first {len(fe['sample_overlapping_ids'])}): "
                  f"{fe['sample_overlapping_ids']}")
    else:
        print("  No violations -- recombination preserved the multitree property.")
    print(f"\n  Total runtime: {runtime_s:.2f}s")


def main():
    parser = argparse.ArgumentParser(description="Per-offspring multitree-violation check (Approach B).")
    parser.add_argument("grg_path", type=Path, help="Path to the .grg file.")
    parser.add_argument("--seed", type=int, default=42,
                        help="numpy random seed for breakpoint generation (default: 42).")
    parser.add_argument("--offspring-per-couple", type=int, default=2,
                        help="Offspring rows produced per mating couple (default: 2).")
    parser.add_argument("--output-json", type=Path, default=None,
                        help="Optional path to dump full per-offspring violation list as JSON.")
    parser.add_argument("--parallel", action="store_true",
                        help="Use ParallelNonDuplicationRecombination (node-aggregated "
                             "discover+commit) instead of the sequential reference.")
    parser.add_argument("--max-workers", type=int, default=None,
                        help="Thread-pool size for --parallel (default: os.cpu_count()).")
    args = parser.parse_args()

    if not args.grg_path.exists():
        print(f"Error: {args.grg_path} not found.", file=sys.stderr)
        sys.exit(1)

    np.random.seed(args.seed)

    t_start = time.perf_counter()

    print(f"Loading {args.grg_path.name}...")
    grg = pygrgl.load_mutable_grg(str(args.grg_path))
    initial_num_nodes = grg.num_nodes
    print(f"  loaded: {initial_num_nodes:,} nodes, {grg.num_edges:,} edges, "
          f"{grg.num_samples:,} samples, {grg.num_mutations:,} mutations")

    engine_label = "parallel" if args.parallel else "sequential"
    print(f"Running one generation of recombination with the {engine_label} engine "
          f"(num_offspring_per_couple={args.offspring_per_couple})...")
    t = time.perf_counter()
    from grg_recombination import simulate_grg_recombination
    if args.parallel:
        from grg_recombination_parallel import ParallelNonDuplicationRecombination
        recomb = ParallelNonDuplicationRecombination(grg, max_workers=args.max_workers)
    else:
        from grg_recombination import NonDuplicationRecombination
        recomb = NonDuplicationRecombination(grg)
    stub = BenchmarkStub(num_offspring_per_couple=args.offspring_per_couple)
    new_offspring_ids, total_bp = simulate_grg_recombination(
        stub, recomb, grg.bp_range, grg.bp_range[1]
    )
    if args.parallel:
        recomb.shutdown()
    print(f"  done in {time.perf_counter() - t:.2f}s "
          f"({len(new_offspring_ids):,} offspring, {total_bp:,} breakpoints, "
          f"{grg.num_nodes - initial_num_nodes:,} new nodes total)")

    print("Computing post-recomb ancestor counts...")
    t = time.perf_counter()
    anc_count = compute_post_recomb_anc_counts(grg)
    print(f"  done in {time.perf_counter() - t:.2f}s")

    print("Running multitree cardinality check on offspring...")
    t = time.perf_counter()
    result = check_offspring(grg, new_offspring_ids, anc_count)
    print(f"  done in {time.perf_counter() - t:.2f}s")

    print_report(args.grg_path, initial_num_nodes, len(new_offspring_ids), result,
                 time.perf_counter() - t_start)

    if args.output_json is not None:
        payload = {
            'grg_path': str(args.grg_path),
            'seed': args.seed,
            'offspring_per_couple': args.offspring_per_couple,
            'initial_num_nodes': initial_num_nodes,
            'num_offspring': len(new_offspring_ids),
            'total_breakpoints': int(total_bp),
            'summary': result['summary'],
            'first_example': result['first_example'],
            'violations': result['violations'],
        }
        with open(args.output_json, 'w') as f:
            json.dump(payload, f, indent=2)
        print(f"\nWrote full violation list to {args.output_json}")


if __name__ == "__main__":
    main()
