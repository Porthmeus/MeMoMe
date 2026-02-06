#!/usr/bin/env python3
"""
Make a 2-panel ORA dotplot (Manual vs MeMoMe/Auto), matching the dotplot style.

Default inputs correspond to:
  - loop-excluding ORA
  - dependency threshold >= 0.8
as produced by `04_ora_subsystems.py`.
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
    if len(q) > 5:
        q = [q[0], q[len(q) // 2], q[-1]]
    return q


def _plot_panel(*, ax, df: pd.DataFrame, title: str, norm, size_scale: float) -> None:
    y = np.arange(len(df))
    x = df["gene_ratio"].to_numpy(dtype=float)
    counts = df["k_high_in_subsystem"].to_numpy(dtype=float)
    fdr = df["fdr_bh"].to_numpy(dtype=float)

    s = np.clip(counts, 1.0, None) * size_scale

    ax.set_axisbelow(True)
    ax.grid(True, which="major", axis="both", color="#e6e6e6", linewidth=1)
    sc = ax.scatter(
        x,
        y,
        s=s,
        c=fdr,
        cmap="coolwarm_r",
        norm=norm,
        edgecolors="none",
        alpha=0.95,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(df["subsystem"].tolist())
    ax.set_xlabel("GeneRatio")
    ax.set_title(title)
    return sc


def plot_pair(*, manual_csv: Path, auto_csv: Path, out_png: Path, out_pdf: Path, top_n: int) -> None:
    plt, LogNorm = _try_import_matplotlib()
    if plt is None:
        print("matplotlib not available; skipping ORA paired plot.")
        return

    manual = _read_ora(manual_csv).sort_values(["fdr_bh", "p_value", "enrichment"], ascending=[True, True, False])
    auto = _read_ora(auto_csv).sort_values(["fdr_bh", "p_value", "enrichment"], ascending=[True, True, False])

    manual = manual.head(top_n).iloc[::-1].reset_index(drop=True)
    auto = auto.head(top_n).iloc[::-1].reset_index(drop=True)

    # Shared color scaling across panels (FDR), log scale.
    fdr_all = pd.concat([manual["fdr_bh"], auto["fdr_bh"]], ignore_index=True).astype(float)
    fdr_pos = fdr_all[(fdr_all > 0) & np.isfinite(fdr_all)].to_numpy()
    if fdr_pos.size:
        norm = LogNorm(vmin=max(float(fdr_pos.min()), 1e-300), vmax=float(fdr_pos.max()))
    else:
        norm = None

    # Shared size scaling so marker sizes are comparable between panels.
    counts_all = pd.concat([manual["k_high_in_subsystem"], auto["k_high_in_subsystem"]], ignore_index=True).astype(float)
    c_max = float(np.nanmax(np.clip(counts_all.to_numpy(), 1.0, None)))
    max_marker_area = 900.0
    size_scale = max_marker_area / c_max if c_max > 0 else 1.0

    fig, axes = plt.subplots(ncols=2, figsize=(16, max(5, 0.45 * max(len(manual), len(auto)))), sharex=True)
    sc0 = _plot_panel(ax=axes[0], df=manual, title="Manual merge (ORA; ≥0.8)", norm=norm, size_scale=size_scale)
    sc1 = _plot_panel(ax=axes[1], df=auto, title="MeMoMe merge (ORA; ≥0.8)", norm=norm, size_scale=size_scale)

    # One shared colorbar.
    cbar = fig.colorbar(sc1, ax=axes.ravel().tolist(), pad=0.02)
    cbar.set_label("FDR (BH)")

    # One shared size legend.
    legend_vals = _choose_size_legend_values(counts_all.to_numpy())
    if legend_vals:
        handles = [
            axes[1].scatter([], [], s=v * size_scale, c="black", alpha=0.9, edgecolors="none", label=str(v))
            for v in legend_vals
        ]
        axes[1].legend(
            handles=handles,
            title="Count",
            loc="center left",
            bbox_to_anchor=(1.05, 0.2),
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
    root = Path("tests/dat/manually_merged_models/gapseq_recon3D/output/ora_subsystems")
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manual",
        type=Path,
        default=root / "results_exclude_loops" / "manual_ora_ge_0.8.csv",
        help="Manual ORA CSV.",
    )
    parser.add_argument(
        "--auto",
        type=Path,
        default=root / "results_exclude_loops" / "auto_ora_ge_0.8.csv",
        help="Auto ORA CSV.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "plots",
        help="Directory to write the paired plot.",
    )
    parser.add_argument("--top-n", type=int, default=20, help="Number of subsystems per panel.")
    args = parser.parse_args()

    out_png = args.output_dir / "manual_vs_auto_ora_ge_0.8_exclude_loops.png"
    out_pdf = args.output_dir / "manual_vs_auto_ora_ge_0.8_exclude_loops.pdf"
    plot_pair(manual_csv=args.manual, auto_csv=args.auto, out_png=out_png, out_pdf=out_pdf, top_n=args.top_n)
    print(f"wrote: {out_png}")


if __name__ == "__main__":
    main()

