#!/usr/bin/env python3
"""
Generate ORA "dot plots" (GeneRatio vs pathway) similar to clusterProfiler outputs.

Inputs:
  ORA result CSVs produced by `04_ora_subsystems.py`, e.g.:
    - ora_subsystems/results_exclude_loops/auto_ora_ge_0.8.csv
    - ora_subsystems/results_include_loops/manual_ora_ge_0.4.csv

Each ORA CSV is expected to contain:
  - subsystem
  - k_high_in_subsystem
  - n_high_total
  - p_value
  - fdr_bh

Outputs:
  For each input ORA CSV, a dot plot is saved as PNG and PDF.

Plot conventions (matching the example style):
  - x axis: GeneRatio = k_high_in_subsystem / n_high_total
  - y axis: subsystem name
  - point size: k_high_in_subsystem ("Count")
  - point color: fdr_bh (BH-adjusted p value), log-scaled (red = small FDR)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def _try_import_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt  # noqa: WPS433
        from matplotlib.colors import LogNorm  # noqa: WPS433

        return plt, LogNorm
    except Exception:  # noqa: BLE001 - best effort plotting
        return None, None


def _savefig(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")


def _read_ora(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    need = {"subsystem", "k_high_in_subsystem", "n_high_total", "p_value", "fdr_bh"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")

    df = df.copy()
    df["subsystem"] = df["subsystem"].astype(str)
    for col in ("k_high_in_subsystem", "n_high_total", "p_value", "fdr_bh"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["k_high_in_subsystem", "n_high_total", "fdr_bh"])
    df = df[df["n_high_total"] > 0]
    df["gene_ratio"] = df["k_high_in_subsystem"] / df["n_high_total"]
    return df


def _choose_size_legend_values(counts: np.ndarray) -> list[int]:
    counts = np.array([c for c in counts if np.isfinite(c) and c > 0], dtype=float)
    if counts.size == 0:
        return []
    q = np.unique(np.round(np.quantile(counts, [0.25, 0.5, 0.75, 1.0])).astype(int))
    q = [int(x) for x in q if x > 0]
    # Limit legend clutter
    if len(q) > 5:
        q = [q[0], q[len(q) // 2], q[-1]]
    return q


def plot_ora_dotplot(*, ora_csv: Path, out_png: Path, out_pdf: Path, top_n: int) -> None:
    plt, LogNorm = _try_import_matplotlib()
    if plt is None:
        print("matplotlib not available; skipping ORA dotplots.")
        return

    df = _read_ora(ora_csv)
    if df.empty:
        print(f"{ora_csv}: empty after filtering; skipping.")
        return

    df = df.sort_values(["fdr_bh", "p_value", "enrichment"], ascending=[True, True, False], na_position="last")
    df = df.head(top_n)

    # Order so most significant is at the top of the plot.
    df = df.iloc[::-1].reset_index(drop=True)

    y = np.arange(len(df))
    x = df["gene_ratio"].to_numpy(dtype=float)
    counts = df["k_high_in_subsystem"].to_numpy(dtype=float)
    fdr = df["fdr_bh"].to_numpy(dtype=float)

    # Size scaling (Count legend uses raw counts; sizes are proportional).
    # Use a dynamic scale so the largest point stays within a consistent visual size,
    # preventing the legend from overflowing when k is large.
    max_marker_area = 900.0  # points^2
    c_max = float(np.nanmax(np.clip(counts, 1.0, None)))
    size_scale = max_marker_area / c_max if c_max > 0 else 1.0
    s = np.clip(counts, 1.0, None) * size_scale

    # Color scaling: use log-normalization so tiny FDR values are distinguishable.
    fdr_pos = fdr[np.isfinite(fdr) & (fdr > 0)]
    if fdr_pos.size == 0:
        norm = None
        vmin = vmax = None
    else:
        vmin = float(np.nanmin(fdr_pos))
        vmax = float(np.nanmax(fdr_pos))
        norm = LogNorm(vmin=max(vmin, 1e-300), vmax=vmax)

    fig, ax = plt.subplots(figsize=(10, max(4, 0.45 * len(df))))
    ax.set_axisbelow(True)
    ax.grid(True, which="major", axis="both", color="#e6e6e6", linewidth=1)

    sc = ax.scatter(
        x,
        y,
        s=s,
        c=fdr,
        cmap="coolwarm_r",  # red=small, blue=large
        norm=norm,
        edgecolors="none",
        alpha=0.95,
    )

    ax.set_yticks(y)
    ax.set_yticklabels(df["subsystem"].tolist())
    ax.set_xlabel("GeneRatio")
    ax.set_ylabel("")

    title = ora_csv.stem.replace("_", " ")
    ax.set_title(title)

    # Colorbar
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("FDR (BH)")

    # Size legend ("Count")
    legend_vals = _choose_size_legend_values(counts)
    if legend_vals:
        handles = [
            ax.scatter([], [], s=v * size_scale, c="black", alpha=0.9, edgecolors="none", label=str(v))
            for v in legend_vals
        ]
        ax.legend(
            handles=handles,
            title="Count",
            loc="center left",
            bbox_to_anchor=(1.20, 0.2),
            frameon=False,
            labelspacing=1.2,
            handletextpad=1.2,
            borderpad=0.6,
            scatterpoints=1,
        )

    _savefig(fig, out_png)
    _savefig(fig, out_pdf)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory containing results_* subfolders with ORA CSVs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "plots",
        help="Directory where dotplots will be written.",
    )
    parser.add_argument("--top-n", type=int, default=20, help="Number of subsystems to plot per ORA CSV.")
    args = parser.parse_args()

    results_dir = args.results_dir
    out_root = args.output_dir

    ora_csvs = sorted(results_dir.glob("results_*/**/*_ora_ge_*.csv"))
    if not ora_csvs:
        raise SystemExit(f"No ORA CSVs found under: {results_dir}")

    for csv_path in ora_csvs:
        rel = csv_path.relative_to(results_dir)
        out_base = (out_root / rel).with_suffix("")
        out_png = out_base.parent / f"{out_base.name}_dotplot.png"
        out_pdf = out_base.parent / f"{out_base.name}_dotplot.pdf"
        plot_ora_dotplot(ora_csv=csv_path, out_png=out_png, out_pdf=out_pdf, top_n=args.top_n)
        print(f"wrote: {out_png}")

    print(f"dotplots written under: {out_root}")


if __name__ == "__main__":
    main()
