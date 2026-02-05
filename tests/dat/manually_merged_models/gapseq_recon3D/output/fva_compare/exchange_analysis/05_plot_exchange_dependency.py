#!/usr/bin/env python3
"""
Plot exchange-reaction dependency distributions (auto vs manual).

Inputs (produced by 04_exchange_dependency.py):
- auto_exchange_dependency.csv (index: canonical_met_id, column: dependency)
- manual_exchange_dependency.csv (index: canonical_met_id, column: dependency)
- auto_manual_exchange_dependency.csv (index: canonical_met_id, columns: dependency_auto/dependency_manual)

Outputs (written to --output-dir):
- exchange_dependency_hist_all.(png|pdf)
- exchange_dependency_hist_overlap.(png|pdf)
- exchange_delta_dependency_hist.(png|pdf)
- exchange_dependency_min_hist_all.(png|pdf)
- exchange_dependency_min_hist_overlap.(png|pdf)
- exchange_delta_dependency_min_hist.(png|pdf)
- exchange_dependency_max_hist_all.(png|pdf)
- exchange_dependency_max_hist_overlap.(png|pdf)
- exchange_delta_dependency_max_hist.(png|pdf)
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

        return plt
    except Exception:  # noqa: BLE001 - best effort plotting
        return None


def _savefig(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")


def _load_dep(path: Path, dep_col: str) -> pd.Series:
    df = pd.read_csv(path, index_col=0)
    if dep_col not in df.columns:
        raise ValueError(f"{path} missing '{dep_col}' column (has: {list(df.columns)})")
    s = df[dep_col].astype(float)
    s.index = s.index.astype(str)
    return s


def _hist_overlay(
    *,
    manual: np.ndarray,
    auto: np.ndarray,
    edges: np.ndarray,
    title: str,
    xlabel: str,
    out_png: Path,
    out_pdf: Path,
) -> None:
    plt = _try_import_matplotlib()
    if plt is None:
        print("matplotlib not available; skipping plot generation.")
        return

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(manual, bins=edges, alpha=0.5, label="Manual", color="#1f77b4", density=True)
    ax.hist(auto, bins=edges, alpha=0.5, label="MeMoMe (Auto)", color="#ff7f0e", density=True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    ax.set_title(title)
    ax.legend(frameon=False)
    _savefig(fig, out_png)
    _savefig(fig, out_pdf)
    plt.close(fig)


def _hist_single(*, values: np.ndarray, edges: np.ndarray, title: str, xlabel: str, out_png: Path, out_pdf: Path) -> None:
    plt = _try_import_matplotlib()
    if plt is None:
        print("matplotlib not available; skipping plot generation.")
        return

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(values, bins=edges, alpha=0.8, color="#444444", density=True)
    ax.axvline(0.0, color="#000000", linewidth=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    ax.set_title(title)
    _savefig(fig, out_png)
    _savefig(fig, out_pdf)
    plt.close(fig)


def _plot_metric(
    *,
    metric_label: str,
    auto_path: Path,
    manual_path: Path,
    intersection_path: Path,
    output_dir: Path,
    bins: int,
    auto_col: str,
    manual_col: str,
    inter_auto_col: str,
    inter_manual_col: str,
    inter_delta_col: str,
    file_prefix: str,
) -> None:
    """Plot histograms for a single dependency metric."""
    auto_dep = _load_dep(auto_path, auto_col).dropna()
    manual_dep = _load_dep(manual_path, manual_col).dropna()

    a_all = auto_dep.to_numpy(dtype=float)
    m_all = manual_dep.to_numpy(dtype=float)

    if len(a_all) and len(m_all):
        x_min = float(np.nanmin([a_all.min(), m_all.min()]))
        x_max = float(np.nanmax([a_all.max(), m_all.max()]))
        edges_all = np.linspace(x_min, x_max, bins + 1)

        _hist_overlay(
            manual=m_all,
            auto=a_all,
            edges=edges_all,
            title=f"Exchange Dependency Distribution ({metric_label}; All)",
            xlabel="Dependency (relative change in EX feasible bound)",
            out_png=output_dir / f"{file_prefix}_hist_all.png",
            out_pdf=output_dir / f"{file_prefix}_hist_all.pdf",
        )

    inter = pd.read_csv(intersection_path, index_col=0)
    inter.index = inter.index.astype(str)
    need = {inter_auto_col, inter_manual_col, inter_delta_col}
    missing = need - set(inter.columns)
    if missing:
        raise ValueError(f"{intersection_path} missing columns: {sorted(missing)}")

    inter = inter.dropna(subset=[inter_auto_col, inter_manual_col, inter_delta_col])
    a_ov = inter[inter_auto_col].astype(float).to_numpy()
    m_ov = inter[inter_manual_col].astype(float).to_numpy()
    d_ov = inter[inter_delta_col].astype(float).to_numpy()

    if len(a_ov) and len(m_ov):
        x_min_ov = float(np.nanmin([a_ov.min(), m_ov.min()]))
        x_max_ov = float(np.nanmax([a_ov.max(), m_ov.max()]))
        edges_ov = np.linspace(x_min_ov, x_max_ov, bins + 1)

        _hist_overlay(
            manual=m_ov,
            auto=a_ov,
            edges=edges_ov,
            title=f"Exchange Dependency Distribution ({metric_label}; Overlap Only)",
            xlabel="Dependency (relative change in EX feasible bound)",
            out_png=output_dir / f"{file_prefix}_hist_overlap.png",
            out_pdf=output_dir / f"{file_prefix}_hist_overlap.pdf",
        )

    if len(d_ov):
        d_min = float(np.nanmin(d_ov))
        d_max = float(np.nanmax(d_ov))
        edges_d = np.linspace(d_min, d_max, bins + 1)
        _hist_single(
            values=d_ov,
            edges=edges_d,
            title=f"Exchange Dependency Delta ({metric_label}; Auto - Manual; Overlap Only)",
            xlabel="Delta dependency (auto - manual)",
            out_png=output_dir / f"{file_prefix}_delta_hist.png",
            out_pdf=output_dir / f"{file_prefix}_delta_hist.pdf",
        )

    print(f"{metric_label}: all manual={len(m_all)} auto={len(a_all)}; overlap={len(inter)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--auto",
        type=Path,
        default=Path(__file__).resolve().parent / "auto_exchange_dependency.csv",
        help="Auto exchange dependency table (canonical_met_id index).",
    )
    parser.add_argument(
        "--manual",
        type=Path,
        default=Path(__file__).resolve().parent / "manual_exchange_dependency.csv",
        help="Manual exchange dependency table (canonical_met_id index).",
    )
    parser.add_argument(
        "--intersection",
        type=Path,
        default=Path(__file__).resolve().parent / "auto_manual_exchange_dependency.csv",
        help="Auto/manual intersection dependency table (canonical_met_id index).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "results",
        help="Directory to write plots.",
    )
    parser.add_argument("--bins", type=int, default=60, help="Histogram bins.")
    args = parser.parse_args()

    outdir = args.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    # Range-based dependency (current default)
    _plot_metric(
        metric_label="range",
        auto_path=args.auto,
        manual_path=args.manual,
        intersection_path=args.intersection,
        output_dir=outdir,
        bins=args.bins,
        auto_col="dependency",
        manual_col="dependency",
        inter_auto_col="dependency_auto",
        inter_manual_col="dependency_manual",
        inter_delta_col="delta_dependency",
        file_prefix="exchange_dependency",
    )

    # Dependency on feasible min (lb)
    _plot_metric(
        metric_label="min",
        auto_path=args.auto,
        manual_path=args.manual,
        intersection_path=args.intersection,
        output_dir=outdir,
        bins=args.bins,
        auto_col="dependency_min",
        manual_col="dependency_min",
        inter_auto_col="dependency_min_auto",
        inter_manual_col="dependency_min_manual",
        inter_delta_col="delta_dependency_min",
        file_prefix="exchange_dependency_min",
    )

    # Dependency on feasible max (ub)
    _plot_metric(
        metric_label="max",
        auto_path=args.auto,
        manual_path=args.manual,
        intersection_path=args.intersection,
        output_dir=outdir,
        bins=args.bins,
        auto_col="dependency_max",
        manual_col="dependency_max",
        inter_auto_col="dependency_max_auto",
        inter_manual_col="dependency_max_manual",
        inter_delta_col="delta_dependency_max",
        file_prefix="exchange_dependency_max",
    )

    print(f"plots written to: {outdir}")


if __name__ == "__main__":
    main()
