"""
haplotype_oracle.py

Per-offspring mutation-correctness verification for GRG recombination.

For each offspring, verifies that the mutations it inherits (by walking up
the GRG from the offspring node to roots) match exactly the mutations
expected from splicing its two parents' haplotypes at the recorded
breakpoints/segments.

The expected mutation set for an offspring with segments
    [(parent_X, end_1), (parent_Y, end_2), ...]
covering [0, end_1), [end_1, end_2), ... is:

    union over segments of {
        m in ancestral_mutations(source_parent)
        where m.position in [seg_start, seg_end)
    }

Run standalone:
    python haplotype_oracle.py path/to/file.grg [--seed 42] \\
            [--offspring-per-couple 2] [--num-generations 1] \\
            [--output-json out.json]
"""

import argparse
import json
import sys
import time
from collections import deque
from pathlib import Path

import numpy as np


def collect_ancestral_mutation_positions(grg, node_id):
    """All mutation positions inherited by node_id (self + all ancestors).

    BFS walk up via get_up_edges. Returns a set of base-pair positions.
    Positions are stable across sort_mutations() (unlike MutationIds).
    """
    visited = set()
    queue = deque([node_id])
    visited.add(node_id)
    positions = set()

    while queue:
        n = queue.popleft()
        for mut_id in grg.get_mutations_for_node(n):
            positions.add(grg.get_mutation_by_id(mut_id).position)
        for parent in grg.get_up_edges(n):
            if parent not in visited:
                visited.add(parent)
                queue.append(parent)

    return positions


def compute_expected_positions(grg, parent_a, parent_b, segments,
                               parent_pos_cache=None):
    """Compute expected mutation positions for an offspring given segments.

    segments: list of (parent_id, end_coord) as produced by
        recombination_intervals. Each segment covers [start, end) of the
        genome where start is the previous segment's end (or 0 for the
        first segment).

    parent_pos_cache: optional dict[int, set[int]] to cache parent
        position sets across siblings sharing the same parents.
    """
    if parent_pos_cache is None:
        parent_pos_cache = {}

    for p in (parent_a, parent_b):
        if p not in parent_pos_cache:
            parent_pos_cache[p] = collect_ancestral_mutation_positions(grg, p)

    pos_a = parent_pos_cache[parent_a]
    pos_b = parent_pos_cache[parent_b]
    parent_positions = {parent_a: pos_a, parent_b: pos_b}

    expected = set()
    start = 0
    for source_parent, end in segments:
        source_pos = parent_positions[source_parent]
        for p in source_pos:
            if start <= p < end:
                expected.add(p)
        start = end

    return expected


def check_offspring_haplotypes(grg, records, verbose=False,
                               max_failure_detail=10):
    """Verify that each offspring carries exactly the expected mutations.

    Args:
        grg: pygrgl.MutableGRG (post-recombination).
        records: list of OffspringRecord (from simulate_grg_recombination).
        verbose: print per-offspring pass/fail.
        max_failure_detail: cap on detailed failure records returned.

    Returns dict with:
        total_checked, num_pass, num_fail, all_pass, failures.
    """
    num_pass = 0
    num_fail = 0
    failures = []
    parent_pos_cache = {}

    for rec in records:
        expected = compute_expected_positions(
            grg, rec.parent_a, rec.parent_b, rec.segments,
            parent_pos_cache=parent_pos_cache,
        )
        actual = collect_ancestral_mutation_positions(grg, rec.offspring_id)

        if expected == actual:
            num_pass += 1
            if verbose:
                print(f"    PASS offspring {rec.offspring_id} "
                      f"(gen {rec.generation})")
        else:
            num_fail += 1
            missing = expected - actual
            extra = actual - expected
            if len(failures) < max_failure_detail:
                failures.append({
                    'offspring_id': int(rec.offspring_id),
                    'parent_a': int(rec.parent_a),
                    'parent_b': int(rec.parent_b),
                    'generation': rec.generation,
                    'num_segments': len(rec.segments),
                    'missing_count': len(missing),
                    'extra_count': len(extra),
                    'missing_positions': sorted(missing)[:20],
                    'extra_positions': sorted(extra)[:20],
                })
            if verbose:
                print(f"    FAIL offspring {rec.offspring_id} "
                      f"(gen {rec.generation}): "
                      f"{len(missing)} missing, {len(extra)} extra")

    return {
        'total_checked': num_pass + num_fail,
        'num_pass': num_pass,
        'num_fail': num_fail,
        'all_pass': num_fail == 0,
        'failures': failures,
    }


def print_report(grg_path, initial_num_nodes, result, runtime_s):
    print(f"\nHaplotype oracle check on {grg_path.name}:")
    print(f"  Initial nodes: {initial_num_nodes:,}")
    print(f"  Offspring checked: {result['total_checked']:,}")
    print(f"  Pass: {result['num_pass']:,}")
    print(f"  Fail: {result['num_fail']:,}")

    if result['failures']:
        f = result['failures'][0]
        print(f"\n  First failure: offspring {f['offspring_id']}")
        print(f"    parents: {f['parent_a']}, {f['parent_b']}")
        print(f"    {f['num_segments']} segments")
        print(f"    {f['missing_count']} missing, {f['extra_count']} extra")
        if f['missing_positions']:
            print(f"    missing positions (first 20): {f['missing_positions']}")
        if f['extra_positions']:
            print(f"    extra positions (first 20): {f['extra_positions']}")
    else:
        print("  All offspring carry the correct mutations.")
    print(f"\n  Total runtime: {runtime_s:.2f}s")


class _BenchmarkStub:
    """Minimal stand-in for breakpoint generation, mirroring
    multitree_check.BenchmarkStub."""

    def __init__(self, num_offspring_per_couple=2):
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


def main():
    parser = argparse.ArgumentParser(
        description="Per-offspring haplotype correctness check.")
    parser.add_argument("grg_path", type=Path, help="Path to the .grg file.")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--offspring-per-couple", type=int, default=2)
    parser.add_argument("--num-generations", type=int, default=1)
    parser.add_argument("--output-json", type=Path, default=None)
    args = parser.parse_args()

    if not args.grg_path.exists():
        print(f"Error: {args.grg_path} not found.", file=sys.stderr)
        sys.exit(1)

    np.random.seed(args.seed)
    import pygrgl
    from grg_recombination import (
        NonDuplicationRecombination,
        OffspringRecord,
        simulate_grg_recombination,
    )

    t_start = time.perf_counter()

    print(f"Loading {args.grg_path.name}...")
    grg = pygrgl.load_mutable_grg(str(args.grg_path))
    initial_num_nodes = grg.num_nodes
    print(f"  loaded: {initial_num_nodes:,} nodes, {grg.num_edges:,} edges, "
          f"{grg.num_samples:,} samples, {grg.num_mutations:,} mutations")

    offspring_ledger = []

    for gen in range(args.num_generations):
        print(f"Running generation {gen + 1} "
              f"(offspring_per_couple={args.offspring_per_couple})...")
        t = time.perf_counter()
        recomb = NonDuplicationRecombination(grg) if gen == 0 else recomb
        stub = _BenchmarkStub(
            num_offspring_per_couple=args.offspring_per_couple)
        offspring_ids, total_bp = simulate_grg_recombination(
            stub, recomb, grg.bp_range, grg.bp_range[1],
            generation_index=gen, offspring_ledger=offspring_ledger,
        )
        print(f"  done in {time.perf_counter() - t:.2f}s "
              f"({len(offspring_ids):,} offspring, {total_bp:,} breakpoints)")

    print(f"\nRunning haplotype oracle check on "
          f"{len(offspring_ledger):,} offspring...")
    t = time.perf_counter()
    result = check_offspring_haplotypes(grg, offspring_ledger, verbose=True)
    oracle_elapsed = time.perf_counter() - t

    print_report(args.grg_path, initial_num_nodes, result,
                 time.perf_counter() - t_start)

    if args.output_json is not None:
        payload = {
            'grg_path': str(args.grg_path),
            'seed': args.seed,
            'offspring_per_couple': args.offspring_per_couple,
            'num_generations': args.num_generations,
            'initial_num_nodes': initial_num_nodes,
            'total_offspring': len(offspring_ledger),
            'oracle_wallclock_s': oracle_elapsed,
            **result,
        }
        with open(args.output_json, 'w') as f:
            json.dump(payload, f, indent=2)
        print(f"\nWrote results to {args.output_json}")


if __name__ == "__main__":
    main()
