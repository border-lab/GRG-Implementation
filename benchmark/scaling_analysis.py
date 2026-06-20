#!/usr/bin/env python3
"""
scaling_analysis.py — Analyze scaling sweep results and project NumPy times.

Reads benchmark CSVs from a scaling sweep directory, fits scaling models
to measured NumPy recombination times, validates via leave-one-out on the
largest measured point, and projects to sizes where NumPy could not run.

Usage:
    python scaling_analysis.py --input-dir ./output/scaling_sweep
    python scaling_analysis.py --input-dir ./output/scaling_sweep --output-dir ./output/scaling_sweep
"""

import argparse
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


ALL_INDIVIDUALS = np.array([4000, 6000, 8000, 12000, 16000, 32000, 64000], dtype=float)
ALL_SNPS = np.array([500_000, 1_000_000], dtype=float)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _parse_target_snps(filename):
    """Parse the intended SNP target from a filename like '4000inds_500k_snps.grg'.

    The actual num_snps in the CSV is a Poisson realization (e.g. 499832),
    not the exact target. This function recovers the round target from the
    filename so we can group rows by intended sweep category.
    """
    m = re.search(r"(\d+)([kKmM])?_snps", filename)
    if not m:
        return None
    val = int(m.group(1))
    suffix = (m.group(2) or "").lower()
    if suffix == "k":
        return val * 1_000
    if suffix == "m":
        return val * 1_000_000
    return val


def load_sweep_results(input_dir):
    """Load all benchmark_recombination_results_*.csv files from input_dir."""
    input_dir = Path(input_dir)
    csv_files = sorted(input_dir.glob("benchmark_recombination_results_*.csv"))

    if not csv_files:
        print(f"No benchmark CSV files found in {input_dir}")
        sys.exit(1)

    int_fields = (
        "num_samples_initial", "num_snps", "num_runs", "num_offspring",
        "nodes_added", "edges_added",
    )
    float_fields = (
        "mean_time_ms", "std_time_ms", "min_time_ms", "max_time_ms",
        "memory_mb", "num_bp", "mean_bp",
        "mean_dead_nodes", "mean_dead_pct", "mean_dead_mutations",
        "mean_dead_edges", "mean_dead_edges_pct",
    )

    rows = []
    for csv_file in csv_files:
        with open(csv_file, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                for key in int_fields:
                    if key in row and row[key]:
                        row[key] = int(float(row[key]))
                for key in float_fields:
                    if key in row and row[key]:
                        row[key] = float(row[key])
                row["num_individuals"] = row["num_samples_initial"] // 2
                # Parse target SNP count from filename (round number used
                # for grouping) instead of actual num_snps (Poisson count).
                target = _parse_target_snps(row.get("file", ""))
                row["snps_target"] = target if target is not None else row["num_snps"]
                rows.append(row)

    print(f"Loaded {len(rows)} rows from {len(csv_files)} CSV files")
    return rows


def separate_by_impl(rows):
    """Split rows into GRG and NumPy lists."""
    grg = [r for r in rows if r["implementation"] == "GRG Native"]
    numpy = [r for r in rows if r["implementation"] == "NumPy Baseline"]
    return grg, numpy


def filter_by_snps(rows, snps):
    """Return rows matching a target SNP category (parsed from filename)."""
    return [r for r in rows if r["snps_target"] == snps]


def extract_xy(rows, x_key="num_individuals", y_key="mean_time_ms"):
    """Extract sorted (x, y, y_err) arrays from rows."""
    if not rows:
        return np.array([]), np.array([]), np.array([])
    pairs = sorted(rows, key=lambda r: r[x_key])
    x = np.array([r[x_key] for r in pairs], dtype=float)
    y = np.array([r[y_key] for r in pairs], dtype=float)
    y_err = np.array([r.get("std_time_ms", 0.0) for r in pairs], dtype=float)
    return x, y, y_err


# ---------------------------------------------------------------------------
# Model fitting
# ---------------------------------------------------------------------------

def _r_squared(y_actual, y_predicted):
    ss_res = np.sum((y_actual - y_predicted) ** 2)
    ss_tot = np.sum((y_actual - np.mean(y_actual)) ** 2)
    if ss_tot == 0:
        return 1.0 if ss_res == 0 else 0.0
    return 1.0 - ss_res / ss_tot


def _fit_linear(x, y):
    coeffs = np.polyfit(x, y, 1)
    pred_fn = lambda xn: np.polyval(coeffs, xn)
    r2 = _r_squared(y, pred_fn(x))
    return {
        "name": "linear",
        "coeffs": {"a": float(coeffs[0]), "b": float(coeffs[1])},
        "formula": f"t = {coeffs[0]:.6e} * x + {coeffs[1]:.4f}",
        "r_squared": float(r2),
        "predict": pred_fn,
    }


def _fit_quadratic(x, y):
    coeffs = np.polyfit(x, y, 2)
    pred_fn = lambda xn: np.polyval(coeffs, xn)
    r2 = _r_squared(y, pred_fn(x))
    return {
        "name": "quadratic",
        "coeffs": {"a": float(coeffs[0]), "b": float(coeffs[1]), "c": float(coeffs[2])},
        "formula": f"t = {coeffs[0]:.6e} * x^2 + {coeffs[1]:.6e} * x + {coeffs[2]:.4f}",
        "r_squared": float(r2),
        "predict": pred_fn,
    }


def _fit_power_law(x, y):
    log_coeffs = np.polyfit(np.log(x), np.log(y), 1)
    k = float(log_coeffs[0])
    a = float(np.exp(log_coeffs[1]))
    pred_fn = lambda xn: a * np.power(xn, k)
    r2 = _r_squared(y, pred_fn(x))
    return {
        "name": "power_law",
        "coeffs": {"a": a, "k": k},
        "formula": f"t = {a:.6e} * x^{k:.4f}",
        "r_squared": float(r2),
        "predict": pred_fn,
    }


def fit_all_models(x, y):
    """Fit linear, quadratic (if enough points), and power-law models."""
    models = {}
    if len(x) >= 2:
        models["linear"] = _fit_linear(x, y)
    if len(x) >= 3:
        models["quadratic"] = _fit_quadratic(x, y)
    if len(x) >= 2 and np.all(x > 0) and np.all(y > 0):
        models["power_law"] = _fit_power_law(x, y)
    return models


def best_model_name(models, loo_results=None):
    """Return the best model name.

    When LOO results are available, pick the model with the lowest LOO
    relative error (generalization ability). Fall back to highest R-squared
    when LOO is unavailable or has errors.
    """
    if not models:
        return None
    if loo_results and "error" not in loo_results:
        candidates = {k: v for k, v in loo_results.items() if k in models}
        if candidates:
            return min(candidates, key=lambda k: candidates[k]["rel_error_pct"])
    return max(models, key=lambda k: models[k]["r_squared"])


# ---------------------------------------------------------------------------
# Combined 2D model fitting (all NumPy points together)
# ---------------------------------------------------------------------------

def fit_combined_2d(numpy_rows):
    """Fit t = a*individuals + b*mutations + c using all NumPy data points.

    Returns a dict with coefficients, predict function, R², residuals,
    or None if there are fewer than 3 points.
    """
    if len(numpy_rows) < 3:
        return None

    n = np.array([r["num_individuals"] for r in numpy_rows], dtype=float)
    m = np.array([r["snps_target"] for r in numpy_rows], dtype=float)
    t = np.array([r["mean_time_ms"] for r in numpy_rows], dtype=float)

    X = np.column_stack([n, m, np.ones(len(n))])
    coeffs, residuals, rank, sv = np.linalg.lstsq(X, t, rcond=None)
    a, b, c = float(coeffs[0]), float(coeffs[1]), float(coeffs[2])

    def pred_fn(n_val, m_val):
        n_val = np.asarray(n_val, dtype=float)
        m_val = np.asarray(m_val, dtype=float)
        return a * n_val + b * m_val + c

    t_pred = pred_fn(n, m)
    r2 = float(_r_squared(t, t_pred))

    per_point = []
    for i, row in enumerate(numpy_rows):
        per_point.append({
            "num_individuals": int(n[i]),
            "num_snps": int(m[i]),
            "actual_ms": float(t[i]),
            "predicted_ms": float(t_pred[i]),
            "residual_ms": float(t[i] - t_pred[i]),
        })

    return {
        "name": "combined_2d_linear",
        "coeffs": {"a_per_individual": a, "b_per_mutation": b, "c_constant": c},
        "formula": f"t = {a:.6e} * n + {b:.6e} * m + {c:.4f}",
        "r_squared": r2,
        "predict": pred_fn,
        "per_point_residuals": per_point,
    }


# ---------------------------------------------------------------------------
# Leave-one-out validation
# ---------------------------------------------------------------------------

def leave_one_out(x, y):
    """Drop the largest-x point, fit on the rest, predict it, report error."""
    if len(x) < 3:
        return {"error": "Need >= 3 measured points for LOO validation"}

    order = np.argsort(x)
    x_sorted, y_sorted = x[order], y[order]

    x_train, y_train = x_sorted[:-1], y_sorted[:-1]
    x_held, y_held = float(x_sorted[-1]), float(y_sorted[-1])

    models = fit_all_models(x_train, y_train)
    results = {}
    for name, m in models.items():
        y_pred = float(m["predict"](np.array([x_held]))[0])
        abs_err = abs(y_pred - y_held)
        rel_err = abs_err / y_held * 100 if y_held != 0 else float("inf")
        results[name] = {
            "x_held_out": x_held,
            "y_actual": y_held,
            "y_predicted": y_pred,
            "abs_error_ms": abs_err,
            "rel_error_pct": rel_err,
            "r_squared_on_training": m["r_squared"],
        }

    return results


# ---------------------------------------------------------------------------
# Projection
# ---------------------------------------------------------------------------

def project_numpy_times(models, best_name, x_measured, x_all):
    """Predict times for x values in x_all that have no measurement."""
    x_missing = np.array([v for v in x_all if v not in set(x_measured)])
    if len(x_missing) == 0 or best_name not in models:
        return []

    pred_fn = models[best_name]["predict"]
    y_proj = pred_fn(x_missing)
    return [
        {"num_individuals": float(xv), "projected_time_ms": float(yv)}
        for xv, yv in zip(x_missing, y_proj)
    ]


def numpy_memory_analytical(num_individuals, num_snps):
    """Analytical NumPy matrix memory: mutations * 4n bytes (post-recomb)."""
    return num_snps * 4.0 * num_individuals / (1024 ** 2)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _snps_label(snps):
    if snps >= 1_000_000 and snps % 1_000_000 == 0:
        return f"{int(snps // 1_000_000)}m"
    if snps >= 1_000 and snps % 1_000 == 0:
        return f"{int(snps // 1_000)}k"
    return str(int(snps))


def _inds_ticks(x_vals):
    """Human-readable tick labels for individual counts."""
    labels = []
    for v in x_vals:
        if v >= 1000:
            labels.append(f"{int(v // 1000)}k")
        else:
            labels.append(str(int(v)))
    return labels


def plot_time_vs_individuals(
    grg_rows, numpy_rows, projections, models, best_name, snps,
    output_path, loo_results=None,
):
    """Figure: time vs individuals at a fixed SNP count."""
    if not HAS_MATPLOTLIB:
        return

    fig, ax = plt.subplots(figsize=(10, 6))
    snps_lbl = _snps_label(snps)

    # GRG measured
    grg_at = filter_by_snps(grg_rows, snps)
    gx, gy, gy_err = extract_xy(grg_at)
    if len(gx):
        ax.errorbar(gx, gy, yerr=gy_err, marker="o", capsize=4,
                     label="GRG Native (measured)", color="#2196F3", linewidth=2)

    # NumPy measured
    np_at = filter_by_snps(numpy_rows, snps)
    nx, ny, ny_err = extract_xy(np_at)
    if len(nx):
        ax.errorbar(nx, ny, yerr=ny_err, marker="s", capsize=4,
                     label="NumPy Baseline (measured)", color="#FF9800", linewidth=2)

    # NumPy projected
    proj_at = [p for p in projections if p.get("num_snps") == snps]
    if proj_at:
        px = np.array([p["num_individuals"] for p in proj_at])
        py = np.array([p["projected_time_ms"] for p in proj_at])
        if len(nx):
            bridge_x = np.array([nx[-1], px[0]])
            bridge_y = np.array([ny[-1], py[0]])
            ax.plot(bridge_x, bridge_y, "--", color="#FF9800", alpha=0.5, linewidth=1.5)
        ax.plot(px, py, "s", markersize=8, markerfacecolor="none",
                markeredgecolor="#FF9800", markeredgewidth=2)
        ax.plot(px, py, "--", color="#FF9800", alpha=0.5, linewidth=1.5,
                label="NumPy Baseline (projected)")

    # Fitted curve (faint, full range)
    if best_name and best_name in models:
        x_curve = np.linspace(ALL_INDIVIDUALS[0], ALL_INDIVIDUALS[-1], 200)
        y_curve = models[best_name]["predict"](x_curve)
        ax.plot(x_curve, y_curve, ":", color="#FF9800", alpha=0.3, linewidth=1,
                label=f"Fit: {models[best_name]['formula']}")

    # LOO annotation
    if loo_results and best_name in loo_results:
        loo = loo_results[best_name]
        ax.annotate(
            f"LOO: {loo['rel_error_pct']:.1f}% error",
            xy=(loo["x_held_out"], loo["y_actual"]),
            xytext=(15, 15), textcoords="offset points",
            fontsize=9, color="gray",
            arrowprops=dict(arrowstyle="->", color="gray", lw=0.8),
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of Individuals", fontsize=12)
    ax.set_ylabel("Mean Recombination Time (ms)", fontsize=12)
    ax.set_title(f"Recombination Time vs Individuals ({snps_lbl} SNPs)", fontsize=14)
    ax.set_xticks(ALL_INDIVIDUALS)
    ax.set_xticklabels(_inds_ticks(ALL_INDIVIDUALS))
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"  Saved {output_path}")


def plot_time_vs_snps(grg_rows, numpy_rows, output_path):
    """Figure: time vs SNPs, one line per individual count."""
    if not HAS_MATPLOTLIB:
        return

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(ALL_INDIVIDUALS)))

    for i, n_ind in enumerate(ALL_INDIVIDUALS):
        ind_lbl = _inds_ticks([n_ind])[0]

        grg_at = [r for r in grg_rows if r["num_individuals"] == n_ind]
        gx, gy, gy_err = extract_xy(grg_at, x_key="snps_target")
        if len(gx):
            ax.errorbar(gx, gy, yerr=gy_err, marker="o", capsize=4,
                         color=colors[i], linewidth=2,
                         label=f"GRG {ind_lbl} inds")

        np_at = [r for r in numpy_rows if r["num_individuals"] == n_ind]
        nx, ny, ny_err = extract_xy(np_at, x_key="snps_target")
        if len(nx):
            ax.errorbar(nx, ny, yerr=ny_err, marker="s", capsize=4,
                         color=colors[i], linewidth=2, linestyle="--",
                         label=f"NumPy {ind_lbl} inds")

    ax.set_yscale("log")
    ax.set_xlabel("Number of SNPs", fontsize=12)
    ax.set_ylabel("Mean Recombination Time (ms)", fontsize=12)
    ax.set_title("Recombination Time vs Mutations", fontsize=14)
    ax.set_xticks(ALL_SNPS)
    ax.set_xticklabels([_snps_label(s) for s in ALL_SNPS])
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=9, ncol=2)
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"  Saved {output_path}")


def plot_memory_vs_individuals(grg_rows, numpy_rows, output_path):
    """Figure: memory vs individuals for both SNP counts."""
    if not HAS_MATPLOTLIB:
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    for snps in ALL_SNPS:
        snps_lbl = _snps_label(snps)

        # GRG measured RSS
        grg_at = filter_by_snps(grg_rows, snps)
        gx, gy_mem, _ = extract_xy(grg_at, y_key="memory_mb")
        if len(gx):
            ax.plot(gx, gy_mem / 1024, marker="o", linewidth=2,
                    label=f"GRG RSS ({snps_lbl} SNPs)")

        # NumPy analytical memory (full range)
        mem_all = np.array([numpy_memory_analytical(n, snps) for n in ALL_INDIVIDUALS])
        ax.plot(ALL_INDIVIDUALS, mem_all / 1024, marker="s", linewidth=2,
                linestyle="--", label=f"NumPy matrix ({snps_lbl} SNPs)")

    # 128 GB system limit
    ax.axhline(y=128, color="red", linestyle=":", linewidth=1.5, alpha=0.7,
               label="System limit (128 GB)")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of Individuals", fontsize=12)
    ax.set_ylabel("Memory (GB)", fontsize=12)
    ax.set_title("Memory Usage vs Individuals", fontsize=14)
    ax.set_xticks(ALL_INDIVIDUALS)
    ax.set_xticklabels(_inds_ticks(ALL_INDIVIDUALS))
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"  Saved {output_path}")


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def save_summary_csv(grg_rows, numpy_rows, projections, output_path):
    """Write a unified CSV with measured + projected rows."""
    fieldnames = [
        "num_individuals", "num_snps", "implementation",
        "mean_time_ms", "std_time_ms", "memory_mb", "source",
    ]

    rows_out = []
    for r in grg_rows:
        rows_out.append({
            "num_individuals": r["num_individuals"],
            "num_snps": int(r["snps_target"]),
            "implementation": "GRG Native",
            "mean_time_ms": f"{r['mean_time_ms']:.4f}",
            "std_time_ms": f"{r.get('std_time_ms', 0.0):.4f}",
            "memory_mb": f"{r['memory_mb']:.1f}",
            "source": "measured",
        })
    for r in numpy_rows:
        rows_out.append({
            "num_individuals": r["num_individuals"],
            "num_snps": int(r["snps_target"]),
            "implementation": "NumPy Baseline",
            "mean_time_ms": f"{r['mean_time_ms']:.4f}",
            "std_time_ms": f"{r.get('std_time_ms', 0.0):.4f}",
            "memory_mb": f"{r['memory_mb']:.1f}",
            "source": "measured",
        })
    for p in projections:
        n_ind = p["num_individuals"]
        snps = p["num_snps"]
        mem = numpy_memory_analytical(n_ind, snps)
        rows_out.append({
            "num_individuals": int(n_ind),
            "num_snps": int(snps),
            "implementation": "NumPy Baseline",
            "mean_time_ms": f"{p['projected_time_ms']:.4f}",
            "std_time_ms": "NaN",
            "memory_mb": f"{mem:.1f}",
            "source": f"projected ({p.get('model_used', 'unknown')})",
        })

    rows_out.sort(key=lambda r: (r["implementation"], int(r["num_snps"]), int(r["num_individuals"])))

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"  Saved {output_path}")


def save_model_fits_json(all_fits, output_path):
    """Save model fit details, LOO results, and projections to JSON."""
    serializable = {}
    for key, val in all_fits.items():
        entry = {}
        for subkey, subval in val.items():
            if subkey == "models":
                entry["models"] = {
                    name: {k: v for k, v in m.items() if k != "predict"}
                    for name, m in subval.items()
                }
            else:
                entry[subkey] = subval
        serializable[key] = entry

    with open(output_path, "w") as f:
        json.dump(serializable, f, indent=2)

    print(f"  Saved {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def analyze_sweep_axis(numpy_rows, snps, label):
    """Run full analysis for one sweep axis (fixed SNP count, vary individuals).

    Returns (models, best_name, loo_results, projections).
    """
    np_at = filter_by_snps(numpy_rows, snps)
    x, y, _ = extract_xy(np_at)

    print(f"\n{'='*60}")
    print(f"Sweep: {label} (fixed {_snps_label(snps)} SNPs)")
    print(f"  Measured NumPy points: {len(x)}")

    if len(x) < 2:
        print("  Not enough points for model fitting (need >= 2)")
        return {}, None, {}, []

    for xi, yi in zip(x, y):
        print(f"    {_inds_ticks([xi])[0]} inds -> {yi:.2f} ms")

    # Fit models
    models = fit_all_models(x, y)

    # LOO validation (run before model selection so LOO informs the choice)
    loo_results = leave_one_out(x, y)

    # Select best model: prefer lowest LOO error when available, else highest R²
    best = best_model_name(models, loo_results)

    print(f"\n  Model fits:")
    for name, m in models.items():
        marker = " <-- best" if name == best else ""
        print(f"    {name:12s}  R²={m['r_squared']:.6f}  {m['formula']}{marker}")

    if "error" not in loo_results:
        print(f"\n  Leave-one-out (held out: {_inds_ticks([x[-1]])[0]} inds):")
        for name, res in loo_results.items():
            marker = " <-- best" if name == best else ""
            print(f"    {name:12s}  predicted={res['y_predicted']:.2f} ms  "
                  f"actual={res['y_actual']:.2f} ms  "
                  f"error={res['rel_error_pct']:.1f}%{marker}")
    else:
        print(f"\n  LOO: {loo_results['error']}")

    # Project
    measured_inds = set(x.tolist())
    projections = project_numpy_times(models, best, x, ALL_INDIVIDUALS)
    for p in projections:
        p["num_snps"] = snps
        p["model_used"] = best
    if projections:
        print(f"\n  Projections (model: {best}):")
        for p in projections:
            print(f"    {_inds_ticks([p['num_individuals']])[0]} inds -> "
                  f"{p['projected_time_ms']:.2f} ms (projected)")
    else:
        print(f"\n  No projections needed (all points measured)")

    return models, best, loo_results, projections


def main():
    parser = argparse.ArgumentParser(
        description="Analyze scaling sweep results and project NumPy times."
    )
    parser.add_argument(
        "--input-dir", type=Path, required=True,
        help="Directory containing benchmark CSV files from the scaling sweep",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=None,
        help="Directory for output files. Defaults to --input-dir.",
    )
    args = parser.parse_args()

    output_dir = args.output_dir or args.input_dir
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load data
    rows = load_sweep_results(args.input_dir)
    grg_rows, numpy_rows = separate_by_impl(rows)
    print(f"  GRG rows: {len(grg_rows)}, NumPy rows: {len(numpy_rows)}")

    # 2. Combined 2D model (all NumPy points together)
    combined = fit_combined_2d(numpy_rows)
    all_projections = []

    if combined:
        print(f"\n{'='*60}")
        print("Combined 2D Model (all NumPy points)")
        print(f"  {combined['formula']}")
        print(f"  R² = {combined['r_squared']:.6f}")
        coeffs = combined["coeffs"]
        print(f"  Per-individual contribution: {coeffs['a_per_individual']:.4f} ms/individual")
        print(f"  Per-mutation contribution:   {coeffs['b_per_mutation']:.6f} ms/mutation")
        print(f"  Constant baseline:           {coeffs['c_constant']:.2f} ms")

        print(f"\n  Per-point residuals:")
        for pt in combined["per_point_residuals"]:
            ind_lbl = _inds_ticks([pt["num_individuals"]])[0]
            snp_lbl = _snps_label(pt["num_snps"])
            print(f"    {ind_lbl}/{snp_lbl}: actual={pt['actual_ms']:.1f}  "
                  f"predicted={pt['predicted_ms']:.1f}  "
                  f"residual={pt['residual_ms']:+.1f} ms")

        # Project to all unmeasured (individuals, mutations) combos
        measured_pairs = {
            (r["num_individuals"], int(r["snps_target"])) for r in numpy_rows
        }
        print(f"\n  Projections:")
        for n_ind in ALL_INDIVIDUALS:
            for snps in ALL_SNPS:
                if (n_ind, int(snps)) in measured_pairs:
                    continue
                t_proj = float(combined["predict"](n_ind, snps))
                all_projections.append({
                    "num_individuals": float(n_ind),
                    "num_snps": float(snps),
                    "projected_time_ms": t_proj,
                    "model_used": "combined_2d_linear",
                })
                ind_lbl = _inds_ticks([n_ind])[0]
                snp_lbl = _snps_label(snps)
                print(f"    {ind_lbl}/{snp_lbl} -> {t_proj:.2f} ms")
    else:
        print("\n  Not enough NumPy points for combined 2D model (need >= 3)")

    # 3. Per-sweep analysis (secondary, for comparison)
    all_fits = {}

    for snps in ALL_SNPS:
        snps_lbl = _snps_label(snps)
        models, best, loo, projections = analyze_sweep_axis(
            numpy_rows, snps, f"vary_individuals_{snps_lbl}",
        )
        all_fits[f"individuals_at_{snps_lbl}"] = {
            "models": models,
            "best_model": best,
            "loo_validation": loo,
            "projections": projections,
        }

    if combined:
        all_fits["combined_2d"] = {
            "name": combined["name"],
            "coeffs": combined["coeffs"],
            "formula": combined["formula"],
            "r_squared": combined["r_squared"],
            "per_point_residuals": combined["per_point_residuals"],
            "projections": [
                {k: v for k, v in p.items() if k != "predict"}
                for p in all_projections
            ],
        }

    # 4. Generate figures
    if not HAS_MATPLOTLIB:
        print("\nSkipping figures (matplotlib not installed)")
    else:
        print(f"\nGenerating figures...")

        for snps in ALL_SNPS:
            snps_lbl = _snps_label(snps)
            proj_at_snps = [p for p in all_projections if p.get("num_snps") == snps]
            combined_models = {}
            combined_best = None
            if combined:
                combined_pred_at_snps = lambda xn, _s=snps, _c=combined: _c["predict"](xn, _s)
                combined_models["combined_2d"] = {
                    "predict": combined_pred_at_snps,
                    "formula": combined["formula"],
                    "r_squared": combined["r_squared"],
                }
                combined_best = "combined_2d"
            plot_time_vs_individuals(
                grg_rows, numpy_rows, proj_at_snps,
                combined_models,
                combined_best,
                snps,
                figures_dir / f"time_vs_individuals_{snps_lbl}.png",
            )

        plot_time_vs_snps(grg_rows, numpy_rows,
                          figures_dir / "time_vs_snps.png")

        plot_memory_vs_individuals(grg_rows, numpy_rows,
                                   figures_dir / "memory_vs_individuals.png")

    # 5. Save outputs
    print(f"\nSaving outputs...")
    save_summary_csv(grg_rows, numpy_rows, all_projections,
                     output_dir / "scaling_summary.csv")
    save_model_fits_json(all_fits, output_dir / "scaling_model_fits.json")

    # 6. Final summary
    print(f"\n{'='*60}")
    print("Summary")
    print(f"{'='*60}")
    print(f"  Measured points:  GRG={len(grg_rows)}, NumPy={len(numpy_rows)}")
    print(f"  Projected points: {len(all_projections)}")
    if combined:
        print(f"  Combined 2D model: {combined['formula']}  "
              f"(R²={combined['r_squared']:.6f})")
    for snps in ALL_SNPS:
        snps_lbl = _snps_label(snps)
        fit_info = all_fits.get(f"individuals_at_{snps_lbl}", {})
        best = fit_info.get("best_model")
        if best and best in fit_info.get("models", {}):
            m = fit_info["models"][best]
            print(f"  Per-sweep at {snps_lbl} SNPs: {best} "
                  f"(R²={m['r_squared']:.6f})")
    if HAS_MATPLOTLIB:
        print(f"  Figures in: {figures_dir}")
    print(f"  Summary CSV: {output_dir / 'scaling_summary.csv'}")
    print(f"  Model fits:  {output_dir / 'scaling_model_fits.json'}")


if __name__ == "__main__":
    main()
