"""
Parity tests: C++ NonDuplicationRecombiner vs Python NonDuplicationRecombination.

For each .grg file passed in via --grg-dir (or the defaults), construct both
backends on a freshly-loaded copy, run one full generation via
simulate_grg_recombination with a fixed PRNG seed on each, and assert that:

  1. The returned offspring NodeID lists are identical.
  2. The full audit dicts match key-for-key, value-for-value.
  3. The post-recombination GRG state matches (num_nodes, num_edges).
  4. Both backends' audit_check identities pass.

Reports wallclock for both, plus the speedup ratio.

Usage:
    cd benchmark/native
    ../../.venv/bin/python test_parity.py
    ../../.venv/bin/python test_parity.py --grg ../../grg_files/16000inds_1m_snps.grg
"""

import argparse
import json
import os
import sys
import time
from datetime import datetime

import numpy as np
import pygrgl


# Make sure we can import from benchmark/ regardless of cwd.
BENCHMARK_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, BENCHMARK_DIR)

# Parity-test outputs live alongside the benchmark output directory but in a
# dedicated subfolder so they don't collide with benchmark_recombination.py's
# CSV/JSON.
PARITY_OUTPUT_DIR = os.path.join(BENCHMARK_DIR, "output", "parity")

from grg_recombination import (
    NonDuplicationRecombination as PyImpl,
    simulate_grg_recombination,
)
from grg_recombination_native import NonDuplicationRecombination as CppImpl


class _FakeBenchmark:
    """Minimal benchmark stub for simulate_grg_recombination."""

    def __init__(self, num_offspring_per_couple=2, expected_crossovers=1.5):
        self.num_offspring_per_couple = num_offspring_per_couple
        self._expected_crossovers = expected_crossovers

    def get_breakpoints(self, bp_range, expected_crossovers=None):
        if expected_crossovers is None:
            expected_crossovers = self._expected_crossovers
        N = bp_range[1]
        num_bp = np.random.poisson(expected_crossovers)
        if num_bp == 0:
            return [], 0
        bps = sorted(np.random.uniform(0, N, num_bp).astype(int).tolist())
        return bps, num_bp


def run_one(grg_path, seed, num_generations, num_offspring_per_couple):
    """Run num_generations on both backends, return diagnostic dicts.

    Each generation runs under the SAME PRNG seed (re-seed before each call)
    on each backend, so breakpoint sampling and parent shuffling are identical
    across implementations.
    """
    bench = _FakeBenchmark(num_offspring_per_couple=num_offspring_per_couple)

    # --- Python reference ---
    g_py = pygrgl.load_mutable_grg(grg_path)
    py_recomb = PyImpl(g_py, instrument=True)
    py_offspring_per_gen = []
    py_total_t = 0.0
    for gen in range(num_generations):
        np.random.seed(seed + gen)
        t0 = time.perf_counter()
        offspring, _ = simulate_grg_recombination(bench, py_recomb, g_py.bp_range, int(g_py.bp_range[1]))
        py_total_t += time.perf_counter() - t0
        py_offspring_per_gen.append(list(offspring))

    # --- C++ implementation ---
    g_cpp = pygrgl.load_mutable_grg(grg_path)
    cpp_recomb = CppImpl(g_cpp, instrument=True)
    cpp_offspring_per_gen = []
    cpp_total_t = 0.0
    for gen in range(num_generations):
        np.random.seed(seed + gen)
        t0 = time.perf_counter()
        offspring, _ = simulate_grg_recombination(bench, cpp_recomb, g_cpp.bp_range, int(g_cpp.bp_range[1]))
        cpp_total_t += time.perf_counter() - t0
        cpp_offspring_per_gen.append(list(offspring))

    return {
        "py": {
            "offspring_per_gen": py_offspring_per_gen,
            "audit": py_recomb.audit,
            "stats": py_recomb.stats,
            "num_nodes": g_py.num_nodes,
            "num_edges": g_py.num_edges,
            "wallclock": py_total_t,
            "audit_check": py_recomb.audit_check(raise_on_fail=False),
        },
        "cpp": {
            "offspring_per_gen": cpp_offspring_per_gen,
            "audit": cpp_recomb.audit,
            "stats": cpp_recomb.stats,
            "num_nodes": g_cpp.num_nodes,
            "num_edges": g_cpp.num_edges,
            "wallclock": cpp_total_t,
            "audit_check": cpp_recomb.audit_check(raise_on_fail=False),
        },
    }


def compare_and_report(grg_path, result):
    py = result["py"]
    cpp = result["cpp"]
    name = os.path.basename(grg_path)

    print(f"\n=== {name} ===")
    print(f"  wallclock     py={py['wallclock']*1000:9.2f}ms  cpp={cpp['wallclock']*1000:9.2f}ms"
          f"  speedup={py['wallclock'] / max(cpp['wallclock'], 1e-9):.2f}x")
    print(f"  num_nodes     py={py['num_nodes']:>9d}      cpp={cpp['num_nodes']:>9d}")
    print(f"  num_edges     py={py['num_edges']:>9d}      cpp={cpp['num_edges']:>9d}")

    failures = []

    # (1) Offspring lists identical per generation
    for g, (po, co) in enumerate(zip(py["offspring_per_gen"], cpp["offspring_per_gen"])):
        if po != co:
            n_mismatch = sum(1 for a, b in zip(po, co) if a != b)
            failures.append(
                f"gen {g}: offspring lists differ ({n_mismatch}/{len(po)} mismatched, "
                f"py len={len(po)} cpp len={len(co)})"
            )

    # (2) Final audit dict matches
    audit_keys = set(py["audit"].keys()) | set(cpp["audit"].keys())
    audit_diffs = [(k, py["audit"].get(k, "MISSING"), cpp["audit"].get(k, "MISSING"))
                   for k in audit_keys if py["audit"].get(k) != cpp["audit"].get(k)]
    if audit_diffs:
        failures.append(f"audit diffs ({len(audit_diffs)}):")
        for k, pv, cv in audit_diffs[:10]:
            failures.append(f"  {k}: py={pv}, cpp={cv}")

    # (3) GRG state matches
    if py["num_nodes"] != cpp["num_nodes"]:
        failures.append(f"num_nodes mismatch py={py['num_nodes']} cpp={cpp['num_nodes']}")
    if py["num_edges"] != cpp["num_edges"]:
        failures.append(f"num_edges mismatch py={py['num_edges']} cpp={cpp['num_edges']}")

    # (4) Audit identities pass on both
    for backend_name, ck in [("py", py["audit_check"]), ("cpp", cpp["audit_check"])]:
        for ident, r in ck.items():
            if not r["pass"]:
                failures.append(f"{backend_name} audit identity FAIL [{ident}]: "
                                f"lhs={r['lhs']} rhs={r['rhs']}")

    # (5) Stats key sets match (values can differ — timing varies)
    py_stats_keys = set(py["stats"].keys())
    cpp_stats_keys = set(cpp["stats"].keys())
    if py_stats_keys != cpp_stats_keys:
        missing = py_stats_keys - cpp_stats_keys
        extra = cpp_stats_keys - py_stats_keys
        failures.append(f"stats key mismatch: missing in cpp={missing}, extra in cpp={extra}")

    if failures:
        print("  PARITY FAILURES:")
        for f in failures:
            print(f"    {f}")
        return False
    print("  PARITY OK")
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--grg",
        action="append",
        default=None,
        help="Specific .grg file(s) to test. Defaults to grg_files/{50inds_1k,500inds_10k} from the repo root.",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--num-generations", type=int, default=1)
    parser.add_argument("--offspring-per-couple", type=int, default=2)
    args = parser.parse_args()

    if args.grg is None:
        repo_root = os.path.abspath(os.path.join(BENCHMARK_DIR, ".."))
        defaults = [
            os.path.join(repo_root, "grg_files", "50inds_1k_snps.grg"),
            os.path.join(repo_root, "grg_files", "500inds_10k_snps.grg"),
        ]
        args.grg = [p for p in defaults if os.path.exists(p)]
        if not args.grg:
            print("No default .grg files found; pass --grg PATH.")
            sys.exit(1)

    os.makedirs(PARITY_OUTPUT_DIR, exist_ok=True)
    run_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    summary_records = []

    overall = True
    for grg_path in args.grg:
        result = run_one(grg_path, args.seed, args.num_generations, args.offspring_per_couple)
        ok = compare_and_report(grg_path, result)
        overall = overall and ok

        # Stash a JSON record per file for later analysis. Keep the structure
        # simple: identifying info, wallclocks, summary stats, audit dicts.
        summary_records.append({
            "grg_file": os.path.basename(grg_path),
            "seed": args.seed,
            "num_generations": args.num_generations,
            "offspring_per_couple": args.offspring_per_couple,
            "parity_passed": ok,
            "py": {
                "wallclock_s": result["py"]["wallclock"],
                "num_nodes": result["py"]["num_nodes"],
                "num_edges": result["py"]["num_edges"],
                "audit": result["py"]["audit"],
            },
            "cpp": {
                "wallclock_s": result["cpp"]["wallclock"],
                "num_nodes": result["cpp"]["num_nodes"],
                "num_edges": result["cpp"]["num_edges"],
                "audit": result["cpp"]["audit"],
            },
            "speedup": result["py"]["wallclock"] / max(result["cpp"]["wallclock"], 1e-9),
        })

    output_path = os.path.join(PARITY_OUTPUT_DIR, f"parity_results_{run_timestamp}.json")
    with open(output_path, "w") as f:
        json.dump({
            "timestamp": run_timestamp,
            "args": {
                "seed": args.seed,
                "num_generations": args.num_generations,
                "offspring_per_couple": args.offspring_per_couple,
                "grg_files": [os.path.basename(p) for p in args.grg],
            },
            "all_passed": overall,
            "records": summary_records,
        }, f, indent=2)
    print(f"\nParity results saved to: {output_path}")

    print()
    print("ALL PARITY TESTS PASSED" if overall else "PARITY FAILURES — see above")
    sys.exit(0 if overall else 1)


if __name__ == "__main__":
    main()
