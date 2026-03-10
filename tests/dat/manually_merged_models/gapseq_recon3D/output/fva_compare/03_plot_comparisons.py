#!/usr/bin/env python3
"""
Plot auto vs manual comparisons (canonicalized, intersection-only).

Inputs:
- auto_manual_host_fva_comparison.csv produced by 02_postprocess_fva_analysis.py

Outputs (written to --output-dir):
- relative_change_hist.(png|pdf)
- relative_change_scatter.(png|pdf)
- delta_rel_change_hist.(png|pdf)
- range_scatter.(png|pdf)
- top_abs_delta_rel_change.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def _try_import_matplotlib():
    try:
        import matplotlib

        # Headless-safe
        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt  # noqa: WPS433

        return plt
    except Exception:  # noqa: BLE001 - best effort for plotting
        return None


def _savefig(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")


def _require_columns(df: pd.DataFrame, cols: list[str]) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def load_comparison(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    _require_columns(
        df,
        [
            "baseline_range_auto",
            "closed_range_auto",
            "relative_change_auto",
            "baseline_range_manual",
            "closed_range_manual",
            "relative_change_manual",
        ],
    )
    df.index = df.index.astype(str)
    return df


def add_delta_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["delta_rel_change"] = out["relative_change_auto"] - out["relative_change_manual"]
    out["abs_delta_rel_change"] = out["delta_rel_change"].abs()
    return out


def plot_relative_change_hist(df: pd.DataFrame, outdir: Path, bins: int) -> None:
    plt = _try_import_matplotlib()
    if plt is None:
        print("matplotlib not available; skipping plot generation.")
        return

    # The scatterplot already shows that the two distributions are very similar; we therefore
    # show the MeMoMe distribution alone to avoid redundancy in the main results figure.
    a = df["relative_change_auto"].to_numpy(dtype=float)

    x_min = float(np.nanmin(a))
    x_max = float(np.nanmax(a))
    # Keep a stable axis for dependency values (expected ~[0, 1]).
    x_min = min(-0.05, x_min)
    x_max = max(1.05, x_max)
    edges = np.linspace(x_min, x_max, bins + 1)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(a, bins=edges, alpha=0.7, label="MeMoMe (Auto)", color="#ff7f0e", density=True)
    ax.axvline(0.8, color="#cc0000", linewidth=2)
    ax.set_xlabel("Dependency")
    ax.set_ylabel("Density")
    ax.set_title("Host Reaction Dependency Distribution")
    ax.legend(frameon=False)

    _savefig(fig, outdir / "relative_change_hist.png")
    _savefig(fig, outdir / "relative_change_hist.pdf")
    plt.close(fig)


def plot_delta_hist(df: pd.DataFrame, outdir: Path, bins: int) -> None:
    plt = _try_import_matplotlib()
    if plt is None:
        print("matplotlib not available; skipping plot generation.")
        return

    d = df["delta_rel_change"].to_numpy(dtype=float)
    x_min = float(np.nanmin(d))
    x_max = float(np.nanmax(d))
    edges = np.linspace(x_min, x_max, bins + 1)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(d, bins=edges, alpha=0.8, color="#444444", density=True)
    ax.axvline(0.0, color="#000000", linewidth=1)
    ax.set_xlabel("Delta dependency (auto - manual)")
    ax.set_ylabel("Density")
    ax.set_title("Dependency Delta Distribution")

    _savefig(fig, outdir / "delta_rel_change_hist.png")
    _savefig(fig, outdir / "delta_rel_change_hist.pdf")
    plt.close(fig)


def plot_dependency_scatter(df: pd.DataFrame, outdir: Path, *, highlight_top_n: int, label_top_n: int) -> None:
    plt = _try_import_matplotlib()
    if plt is None:
        print("matplotlib not available; skipping plot generation.")
        return

    x = df["relative_change_manual"].to_numpy(dtype=float)
    y = df["relative_change_auto"].to_numpy(dtype=float)

    # Smaller figure (same font sizes) for easier inclusion in the thesis.
    fig, ax = plt.subplots(figsize=(3.5, 3.5))
    ax.scatter(x, y, s=8, alpha=0.25, color="#444444", linewidths=0)

    # Highlight top deltas (keep this small for readability)
    top = df.sort_values("abs_delta_rel_change", ascending=False).head(max(highlight_top_n, 1))
    ax.scatter(
        top["relative_change_manual"],
        top["relative_change_auto"],
        s=20,
        alpha=0.9,
        color="#d62728",
        linewidths=0,
    )

    label_n = min(label_top_n, len(top))
    if label_n > 0:
        # Place labels in data coordinates and nudge them to avoid overlaps.
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()

        lim_min = float(np.nanmin([x.min(), y.min()]))
        lim_max = float(np.nanmax([x.max(), y.max()]))
        span = max(lim_max - lim_min, 1e-9)
        dx = 0.01 * span
        dy = 0.01 * span
        step = 0.03 * span

        placed_bboxes = []
        for rid, row in top.head(label_n).iterrows():
            px = float(row["relative_change_manual"])
            py = float(row["relative_change_auto"])

            for attempt in range(50):
                # Prefer placing above the point; if we run out of room, go below.
                if py + dy + attempt * step <= lim_max:
                    ty = py + dy + attempt * step
                else:
                    ty = py - dy - attempt * step

                tx = px + dx
                text = ax.text(tx, ty, rid, fontsize=7, color="#d62728", zorder=5)
                fig.canvas.draw()
                bbox = text.get_window_extent(renderer).expanded(1.05, 1.2)
                if any(bbox.overlaps(other) for other in placed_bboxes):
                    text.remove()
                    continue

                placed_bboxes.append(bbox)
                break

    # y=x
    lim_min = float(np.nanmin([x.min(), y.min()]))
    lim_max = float(np.nanmax([x.max(), y.max()]))
    ax.plot([lim_min, lim_max], [lim_min, lim_max], color="#000000", linewidth=1, alpha=0.6)

    ax.set_xlim(lim_min, lim_max)
    ax.set_ylim(lim_min, lim_max)
    ax.set_xlabel("Manual model dependency")
    ax.set_ylabel("MeMoMe model dependency")
    ax.set_title("Host reactions microbiome-dependency")

    _savefig(fig, outdir / "relative_change_scatter.png")
    _savefig(fig, outdir / "relative_change_scatter.pdf")
    plt.close(fig)


def plot_range_scatter(df: pd.DataFrame, outdir: Path) -> None:
    plt = _try_import_matplotlib()
    if plt is None:
        print("matplotlib not available; skipping plot generation.")
        return

    eps = 1e-6
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharex=True, sharey=True)

    def _panel(ax, xcol: str, ycol: str, title: str) -> None:
        x = np.log10(df[xcol].to_numpy(dtype=float) + eps)
        y = np.log10(df[ycol].to_numpy(dtype=float) + eps)
        ax.scatter(x, y, s=8, alpha=0.25, color="#444444", linewidths=0)
        lim_min = float(np.nanmin([x.min(), y.min()]))
        lim_max = float(np.nanmax([x.max(), y.max()]))
        ax.plot([lim_min, lim_max], [lim_min, lim_max], color="#000000", linewidth=1, alpha=0.6)
        ax.set_title(title)
        ax.set_xlabel("log10(manual range + eps)")

    _panel(axes[0], "baseline_range_manual", "baseline_range_auto", "Baseline ranges")
    _panel(axes[1], "closed_range_manual", "closed_range_auto", "Closed ranges")
    axes[0].set_ylabel("log10(auto range + eps)")

    fig.suptitle("Host Reaction Range Comparison (Overlap Only)")
    _savefig(fig, outdir / "range_scatter.png")
    _savefig(fig, outdir / "range_scatter.pdf")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).resolve().parent / "auto_manual_host_fva_comparison.csv",
        help="Canonicalized intersection comparison CSV (from postprocess).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "results",
        help="Directory to write plots and helper tables.",
    )
    parser.add_argument("--bins", type=int, default=60, help="Histogram bins.")
    parser.add_argument(
        "--highlight-top-n",
        type=int,
        default=10,
        help="Highlight top-N absolute delta reactions on the scatterplot.",
    )
    parser.add_argument(
        "--label-top-n",
        type=int,
        default=10,
        help="Annotate (label) up to top-N highlighted reactions on the scatterplot.",
    )
    parser.add_argument(
        "--top-delta-csv-n",
        type=int,
        default=200,
        help="How many top |delta| reactions to export as CSV.",
    )
    args = parser.parse_args()

    df = load_comparison(args.input)
    df = add_delta_columns(df)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    top = df.sort_values("abs_delta_rel_change", ascending=False).head(max(args.top_delta_csv_n, 1))
    top.to_csv(args.output_dir / "top_abs_delta_rel_change.csv")

    plot_relative_change_hist(df, args.output_dir, bins=args.bins)
    plot_delta_hist(df, args.output_dir, bins=args.bins)
    plot_dependency_scatter(
        df,
        args.output_dir,
        highlight_top_n=args.highlight_top_n,
        label_top_n=args.label_top_n,
    )
    plot_range_scatter(df, args.output_dir)

    # Minimal console summary (useful when running long pipelines)
    top_auto = set(df.sort_values("relative_change_auto", ascending=False).head(30).index)
    top_manual = set(df.sort_values("relative_change_manual", ascending=False).head(30).index)
    print(f"n_reactions={len(df)}")
    print(f"top30_overlap={len(top_auto & top_manual)}")


if __name__ == "__main__":
    main()
