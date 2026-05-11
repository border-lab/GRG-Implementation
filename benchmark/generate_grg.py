#!/usr/bin/env python3
"""
Generate GRG files for benchmarking, either by msprime simulation or by
converting a pre-existing tskit .trees file.

When simulating, number of diploid individuals and (target) number of mutations
are both configurable. Mutation count is controlled either directly via the
per-bp per-generation rate, or by passing a target count which is used to
calibrate the rate from the simulated ancestry's total branch length.

When --input-trees is given, msprime is skipped entirely and the supplied
tree sequence is fed straight into the selected --method conversion path.

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
import tskit


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


def _derive_construct_flags(
    sequence_length: int, segment_size: int, num_jobs: int
) -> tuple[int, int]:
    """Pick `grg construct` -t / -p per DeHaas et al. Methods section
    'Determining Parameters T and P': T = ceil(seq_len / segment_size);
    P bounded by cores (>=) and T (<=)."""
    trees = max(1, round(sequence_length / segment_size))
    parts = min(num_jobs, trees)
    return trees, parts


def generate_grg(
    output_path: Path,
    *,
    num_individuals: int | None = None,
    target_mutations: int | None = None,
    mutation_rate: float | None = None,
    sequence_length: int = 1_000_000,
    Ne: float = 10_000.0,
    recombination_rate: float = 1e-8,
    seed: int | None = None,
    method: str = "ts",
    construct_jobs: int = 1,
    construct_trees: int | None = None,
    construct_parts: int | None = None,
    segment_size: int = 100_000,
    keep_vcf: bool = False,
    input_trees: Path | None = None,
) -> tuple[Path, int]:
    """
    Generate the requested output file. Returns (output_path, mutation_count).

    `method` picks the final stage: "ts" runs `grg convert` on a .trees dump
    (paper's TS-conversion path), "vcf" runs `grg construct` on a .vcf dump
    (paper's BuildShape+MapMutations path; what biobank deployments use), and
    "trees" skips conversion and writes the raw tskit .trees file.

    Source of the tree sequence: if `input_trees` is given, it is loaded via
    tskit and used directly (msprime is skipped, simulation kwargs ignored).
    Otherwise, msprime simulates one — `num_individuals` is required, and
    exactly one of `target_mutations` or `mutation_rate` must be provided.
    """
    if method not in {"ts", "vcf", "trees"}:
        raise ValueError(f"method must be 'ts', 'vcf', or 'trees', got {method!r}")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if input_trees is not None:
        if method == "trees":
            raise ValueError(
                "method='trees' is incompatible with input_trees (would just "
                "copy the input file)."
            )
        ts = tskit.load(str(input_trees))
        actual_muts = ts.num_mutations
        print(
            f"  loaded {input_trees}: {ts.num_samples:,} samples, "
            f"{actual_muts:,} mutations, {ts.num_trees:,} trees"
        )
    else:
        if num_individuals is None:
            raise ValueError("num_individuals is required when input_trees is None")
        if (target_mutations is None) == (mutation_rate is None):
            raise ValueError("provide exactly one of target_mutations or mutation_rate")

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

    # Resolve to absolute since the construct branch runs from a tempdir cwd.
    abs_output = output_path.resolve()

    if method == "trees":
        ts.dump(str(abs_output))
    elif method == "ts":
        grg_bin = _grg_binary()
        with tempfile.TemporaryDirectory() as tmpdir:
            trees_path = os.path.join(tmpdir, "sim.trees")
            ts.dump(trees_path)
            subprocess.check_call([grg_bin, "convert", trees_path, str(abs_output)])
    else:
        grg_bin = _grg_binary()
        with tempfile.TemporaryDirectory() as tmpdir:
            if keep_vcf:
                vcf_path = abs_output.with_suffix(".vcf")
            else:
                vcf_path = Path(tmpdir) / "sim.vcf"
            with open(vcf_path, "w") as fh:
                ts.write_vcf(fh)
            derived_t, derived_p = _derive_construct_flags(
                int(ts.sequence_length), segment_size, construct_jobs,
            )
            trees = construct_trees if construct_trees is not None else derived_t
            parts = construct_parts if construct_parts is not None else derived_p
            # cwd=tmpdir keeps `grg construct`'s per-part .grg shards out of
            # the user's cwd; --force skips the unindexed-VCF stop-warning.
            cmd = [
                grg_bin, "construct", str(vcf_path),
                "-j", str(construct_jobs),
                "-t", str(trees),
                "-p", str(parts),
                "--force",
                "-o", str(abs_output),
            ]
            subprocess.check_call(cmd, cwd=tmpdir)

    print(f"  wrote {output_path}")
    return output_path, actual_muts


def _auto_filename(num_individuals: int, num_mutations: int, ext: str = ".grg") -> str:
    """<N>inds_<M>_snps<ext>, with k/m suffixes when the count is a clean multiple."""
    if num_mutations >= 1_000_000 and num_mutations % 1_000_000 == 0:
        snps = f"{num_mutations // 1_000_000}m"
    elif num_mutations >= 1_000 and num_mutations % 1_000 == 0:
        snps = f"{num_mutations // 1_000}k"
    else:
        snps = str(num_mutations)
    return f"{num_individuals}inds_{snps}_snps{ext}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-n", "--num-individuals", type=int, default=None,
        help="Number of diploid individuals (sample nodes = 2*n). "
             "Required unless --input-trees is given.",
    )

    mut = parser.add_mutually_exclusive_group()
    mut.add_argument(
        "-m", "--target-mutations", type=int,
        help="Target mutation count; rate auto-calibrated (actual count is Poisson)",
    )
    mut.add_argument(
        "--mutation-rate", type=float,
        help="Mutation rate per bp per generation (used directly)",
    )
    parser.add_argument(
        "--input-trees", type=Path, default=None,
        help="Path to an existing tskit .trees file. When set, msprime is "
             "skipped and this tree sequence is fed straight into --method. "
             "Mutually exclusive with -n / -m / --mutation-rate / simulation "
             "knobs; not allowed with --method trees.",
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
    parser.add_argument("--method", choices=["ts", "vcf", "trees"], default="ts",
                        help="Pipeline: 'ts' runs `grg convert` on a .trees dump; "
                             "'vcf' runs `grg construct` on a .vcf dump; "
                             "'trees' writes only the tskit .trees file (no GRG). "
                             "Default: ts")
    parser.add_argument("--construct-jobs", type=int, default=1,
                        help="Threads passed to `grg construct -j` (--method vcf). "
                             "Default: 1")
    parser.add_argument("--construct-trees", type=int, default=None,
                        help="Override for `grg construct -t` (--method vcf). "
                             "Default: derived as seq_len / --segment-size.")
    parser.add_argument("--construct-parts", type=int, default=None,
                        help="Override for `grg construct -p` (--method vcf). "
                             "Default: min(--construct-jobs, derived T).")
    parser.add_argument("--segment-size", type=int, default=100_000,
                        help="Target segment size in bp used to derive "
                             "`grg construct -t` (--method vcf). Paper "
                             "recommends 50_000-150_000. Default: 100_000")
    parser.add_argument("--keep-vcf", action="store_true",
                        help="Keep the intermediate VCF next to the output GRG "
                             "(--method vcf only).")

    args = parser.parse_args()

    if args.input_trees is not None:
        if args.method == "trees":
            parser.error("--method trees is incompatible with --input-trees "
                         "(would just copy the input file).")
        sim_args_set = [
            ("-n/--num-individuals", args.num_individuals),
            ("-m/--target-mutations", args.target_mutations),
            ("--mutation-rate", args.mutation_rate),
        ]
        conflicts = [name for name, val in sim_args_set if val is not None]
        if conflicts:
            parser.error(
                f"--input-trees cannot be combined with simulation args: "
                f"{', '.join(conflicts)}"
            )
    else:
        if args.num_individuals is None:
            parser.error("-n/--num-individuals is required (or pass --input-trees).")
        if args.target_mutations is None and args.mutation_rate is None:
            parser.error(
                "one of -m/--target-mutations or --mutation-rate is required "
                "(or pass --input-trees)."
            )

    out_ext = ".trees" if args.method == "trees" else ".grg"
    if args.output is None:
        if args.input_trees is not None:
            args.output = args.output_dir / f"{args.input_trees.stem}{out_ext}"
        elif args.target_mutations is not None:
            args.output = args.output_dir / _auto_filename(args.num_individuals, args.target_mutations, out_ext)
        else:
            # Actual count is Poisson under --mutation-rate; rename after we know it.
            args.output = args.output_dir / f"{args.num_individuals}inds_pending{out_ext}"

    out_path, actual = generate_grg(
        output_path=args.output,
        num_individuals=args.num_individuals,
        target_mutations=args.target_mutations,
        mutation_rate=args.mutation_rate,
        sequence_length=args.sequence_length,
        Ne=args.Ne,
        recombination_rate=args.recombination_rate,
        seed=args.seed,
        method=args.method,
        construct_jobs=args.construct_jobs,
        construct_trees=args.construct_trees,
        construct_parts=args.construct_parts,
        segment_size=args.segment_size,
        keep_vcf=args.keep_vcf,
        input_trees=args.input_trees,
    )

    if (
        args.input_trees is None
        and args.target_mutations is None
        and out_path.name.endswith(f"_pending{out_ext}")
    ):
        new_path = out_path.parent / _auto_filename(args.num_individuals, actual, out_ext)
        out_path.rename(new_path)
        print(f"  renamed -> {new_path}")


if __name__ == "__main__":
    main()
