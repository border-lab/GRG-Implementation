#!/usr/bin/env python3
"""
Generate GRG files via msprime simulation for benchmarking.

Number of diploid individuals and (target) number of mutations are both
configurable. Mutation count is controlled either directly via the per-bp
per-generation rate, or by passing a target count which is used to calibrate
the rate from the simulated ancestry's total branch length.

Filenames follow the <N>inds_<M>_snps.grg convention recognized by
benchmark_recombination._parse_expected_size, so generated files can be
dropped straight into ./grg_files for the benchmarker to pick up.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import tempfile
from pathlib import Path
from shutil import which

import msprime


def _grg_binary() -> str:
    grg_bin = which("grg")
    if grg_bin is None:
        raise RuntimeError(
            "Could not find 'grg' on PATH. Activate the env that has pygrgl installed."
        )
    return grg_bin


def _total_branch_length(ts) -> float:
    """Integral of branch length across the genome: sum over edges of
    (parent_time - child_time) * span. Multiplied by a mutation rate, this
    gives the expected number of mutations sim_mutations will draw."""
    edges = ts.tables.edges
    times = ts.tables.nodes.time
    spans = edges.right - edges.left
    return float(((times[edges.parent] - times[edges.child]) * spans).sum())


def generate_grg(
    output_path: Path,
    *,
    num_individuals: int,
    target_mutations: int | None = None,
    mutation_rate: float | None = None,
    sequence_length: int = 1_000_000,
    Ne: float = 10_000.0,
    recombination_rate: float = 1e-8,
    seed: int | None = None,
) -> tuple[Path, int]:
    """
    Simulate ancestry+mutations with msprime and write the result as a .grg
    file via the `grg convert` CLI. Returns (output_path, actual_mutation_count).

    Exactly one of `target_mutations` or `mutation_rate` must be provided.
    """
    if (target_mutations is None) == (mutation_rate is None):
        raise ValueError("provide exactly one of target_mutations or mutation_rate")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    ts = msprime.sim_ancestry(
        samples=num_individuals,
        sequence_length=sequence_length,
        population_size=Ne,
        recombination_rate=recombination_rate,
        random_seed=seed,
    )

    if mutation_rate is None:
        L = _total_branch_length(ts)
        if L <= 0:
            raise RuntimeError("Tree sequence has zero total branch length.")
        mutation_rate = target_mutations / L
        print(
            f"  calibrated mutation_rate={mutation_rate:.3e} "
            f"(target {target_mutations:,} muts / total branch {L:.3e})"
        )

    mut_seed = None if seed is None else seed + 1
    ts = msprime.sim_mutations(ts, rate=mutation_rate, random_seed=mut_seed)
    actual_muts = ts.num_mutations
    print(
        f"  msprime: {ts.num_samples:,} samples, "
        f"{actual_muts:,} mutations, {ts.num_trees:,} trees"
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        trees_path = os.path.join(tmpdir, "sim.trees")
        ts.dump(trees_path)
        subprocess.check_call([_grg_binary(), "convert", trees_path, str(output_path)])

    print(f"  wrote {output_path}")
    return output_path, actual_muts


def _auto_filename(num_individuals: int, num_mutations: int) -> str:
    """<N>inds_<M>_snps.grg, with k/m suffixes when the count is a clean multiple."""
    if num_mutations >= 1_000_000 and num_mutations % 1_000_000 == 0:
        snps = f"{num_mutations // 1_000_000}m"
    elif num_mutations >= 1_000 and num_mutations % 1_000 == 0:
        snps = f"{num_mutations // 1_000}k"
    else:
        snps = str(num_mutations)
    return f"{num_individuals}inds_{snps}_snps.grg"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-n", "--num-individuals", type=int, required=True,
        help="Number of diploid individuals (sample nodes = 2*n)",
    )

    mut = parser.add_mutually_exclusive_group(required=True)
    mut.add_argument(
        "-m", "--target-mutations", type=int,
        help="Target mutation count; rate auto-calibrated (actual count is Poisson)",
    )
    mut.add_argument(
        "--mutation-rate", type=float,
        help="Mutation rate per bp per generation (used directly)",
    )

    parser.add_argument("--sequence-length", type=int, default=1_000_000,
                        help="Sequence length in bp (default: 1e6)")
    parser.add_argument("--Ne", type=float, default=10_000.0,
                        help="Effective population size (default: 10000)")
    parser.add_argument("--recombination-rate", type=float, default=1e-8,
                        help="Recombination rate per bp per generation (default: 1e-8)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="Output path. If omitted, auto-named under --output-dir.")
    parser.add_argument("--output-dir", type=Path, default=Path("./grg_files"),
                        help="Directory for auto-named output (default: ./grg_files)")

    args = parser.parse_args()

    if args.output is None:
        if args.target_mutations is not None:
            args.output = args.output_dir / _auto_filename(args.num_individuals, args.target_mutations)
        else:
            # Actual count is Poisson under --mutation-rate; rename after we know it.
            args.output = args.output_dir / f"{args.num_individuals}inds_pending.grg"

    out_path, actual = generate_grg(
        output_path=args.output,
        num_individuals=args.num_individuals,
        target_mutations=args.target_mutations,
        mutation_rate=args.mutation_rate,
        sequence_length=args.sequence_length,
        Ne=args.Ne,
        recombination_rate=args.recombination_rate,
        seed=args.seed,
    )

    if args.target_mutations is None and out_path.name.endswith("_pending.grg"):
        new_path = out_path.parent / _auto_filename(args.num_individuals, actual)
        out_path.rename(new_path)
        print(f"  renamed -> {new_path}")


if __name__ == "__main__":
    main()
