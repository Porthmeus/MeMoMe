#!/usr/bin/env python3
"""Plot FVA comparison results produced by 01_run_fva_analysis.py."""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import pandas as pd


def _strip_reaction_compartment(rxn_id: str) -> str:
    rid = rxn_id
    if rid.startswith("R_"):
        rid = rid[2:]
    if rid.endswith("_t"):
        rid = rid[:-2]
    rid = re.sub(r"\([A-Za-z0-9_]+\)$", "", rid)
    rid = re.sub(r"\[[A-Za-z0-9_]+\]$", "", rid)
    rid = re.sub(r"_[a-z]\d?$", "", rid)
    return rid


def _canonical_reaction_id(rxn_id: str) -> str:
    rid = rxn_id
    if rid.startswith("R_"):
        rid = rid[2:]
    for prefix in ("TR_model1_", "TR_model2_"):
        if rid.startswith(prefix):
            return _strip_reaction_compartment(rid[len(prefix) :])

    for sp in ("TR_", "IEX_"):
        if rid.startswith(sp):
            rest = rid[len(sp) :]
            for mp in ("model1_", "model2_", "H_", "M_"):
                if rest.startswith(mp):
                    return _strip_reaction_compartment(rest[len(mp) :])
            return _strip_reaction_compartment(rest)

    for sp in ("EX_", "DM_", "SK_", "sink_"):
        if rid.startswith(sp):
            rest = rid[len(sp) :]
            for mp in ("model1_", "model2_", "H_", "M_"):
                if rest.startswith(mp):
                    return _strip_reaction_compartment(sp + rest[len(mp) :])
            return _strip_reaction_compartment(rid)

    for prefix in ("model1_", "model2_", "H_", "M_"):
        if rid.startswith(prefix):
            return _strip_reaction_compartment(rid[len(prefix) :])

    return _strip_reaction_compartment(rid)


def _plot_distribution(
    df: pd.DataFrame,
    output_path: Path,
    title: str,
    xlim: tuple[float, float] | None,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot generation.")
        return

    values = df["relative_change"].dropna()
    if values.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(values, bins=10, color="#4C72B0", alpha=0.85)
    ax.set_title(title)
    ax.set_xlabel("Relative change")
    ax.set_ylabel("Count")
    if xlim is not None:
        ax.set_xlim(xlim)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def _plot_top_changes(
    df: pd.DataFrame,
    output_path: Path,
    title: str,
    top_n: int,
    xlim: float | None,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot generation.")
        return

    values = df["relative_change"].dropna().abs()
    if values.empty:
        return
    top = values.sort_values(ascending=False).head(top_n)
    fig, ax = plt.subplots(figsize=(10, max(4, top_n * 0.25)))
    ax.barh(top.index.astype(str), top.values, color="#DD8452")
    ax.set_title(title)
    ax.set_xlabel("Absolute relative change")
    ax.invert_yaxis()
    if xlim is not None:
        ax.set_xlim(0, xlim)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def _plot_top_overlap(
    auto_df: pd.DataFrame,
    manual_df: pd.DataFrame,
    output_path: Path,
    top_n: int,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot generation.")
        return

    auto_top = {
        _canonical_reaction_id(rid)
        for rid in auto_df["relative_change"].abs().dropna().sort_values(ascending=False).head(top_n).index
    }
    manual_top = {
        _canonical_reaction_id(rid)
        for rid in manual_df["relative_change"].abs().dropna().sort_values(ascending=False).head(top_n).index
    }

    overlap = len(auto_top & manual_top)
    auto_only = len(auto_top - manual_top)
    manual_only = len(manual_top - auto_top)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(["Auto only", "Overlap", "Manual only"], [auto_only, overlap, manual_only], color="#55A868")
    ax.set_title("Top-change reaction overlap")
    ax.set_ylabel("Count")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def _plot_overlap_all(
    auto_df: pd.DataFrame,
    manual_df: pd.DataFrame,
    output_path: Path,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot generation.")
        return

    auto_ids = {_canonical_reaction_id(rid) for rid in auto_df.index}
    manual_ids = {_canonical_reaction_id(rid) for rid in manual_df.index}
    overlap = len(auto_ids & manual_ids)
    auto_only = len(auto_ids - manual_ids)
    manual_only = len(manual_ids - auto_ids)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(["Auto only", "Overlap", "Manual only"], [auto_only, overlap, manual_only], color="#4C72B0")
    ax.set_title("Overall reaction overlap")
    ax.set_ylabel("Count")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def _plot_relative_change_scatter(
    auto_df: pd.DataFrame,
    manual_df: pd.DataFrame,
    output_path: Path,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot generation.")
        return

    auto_map = {}
    for rid, value in auto_df["relative_change"].items():
        if pd.isna(value):
            continue
        auto_map.setdefault(_canonical_reaction_id(rid), []).append(float(value))

    manual_map = {}
    for rid, value in manual_df["relative_change"].items():
        if pd.isna(value):
            continue
        manual_map.setdefault(_canonical_reaction_id(rid), []).append(float(value))

    overlap = sorted(set(auto_map) & set(manual_map))
    if not overlap:
        return

    xs = [sum(auto_map[rid]) / len(auto_map[rid]) for rid in overlap]
    ys = [sum(manual_map[rid]) / len(manual_map[rid]) for rid in overlap]
    max_abs = max((abs(v) for v in xs + ys), default=1.0)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(xs, ys, alpha=0.5, s=10, color="#55A868", edgecolors="none")
    ax.plot([-max_abs, max_abs], [-max_abs, max_abs], color="#999999", linestyle="--", linewidth=1)
    ax.set_xlim(-max_abs, max_abs)
    ax.set_ylim(-max_abs, max_abs)
    ax.set_xlabel("Auto relative change")
    ax.set_ylabel("Manual relative change")
    ax.set_title("Relative change per reaction (overlap)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot FVA comparison results.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare"),
    )
    parser.add_argument(
        "--plot-top-n",
        type=int,
        default=30,
        help="Number of top-changing reactions to plot.",
    )
    args = parser.parse_args()

    auto_path = args.output_dir / "auto_host_fva_comparison.csv"
    manual_path = args.output_dir / "manual_host_fva_comparison.csv"

    if auto_path.exists():
        auto_df = pd.read_csv(auto_path, index_col=0)
    else:
        auto_df = None

    if manual_path.exists():
        manual_df = pd.read_csv(manual_path, index_col=0)
    else:
        manual_df = None

    combined = []
    for df in (auto_df, manual_df):
        if df is None:
            continue
        values = df["relative_change"].replace([float("inf"), float("-inf")], pd.NA).dropna()
        combined.append(values)
    if combined:
        combined_values = pd.concat(combined)
        max_abs = combined_values.abs().max()
        xlim = (-max_abs, max_abs)
    else:
        xlim = None
        max_abs = None

    if auto_df is not None:
        _plot_distribution(
            auto_df,
            args.output_dir / "auto_relative_change_hist.png",
            "Auto relative change",
            xlim,
        )
        _plot_top_changes(
            auto_df,
            args.output_dir / "auto_top_relative_change.png",
            "Auto top changes",
            args.plot_top_n,
            max_abs,
        )

    if manual_df is not None:
        _plot_distribution(
            manual_df,
            args.output_dir / "manual_relative_change_hist.png",
            "Manual relative change",
            xlim,
        )
        _plot_top_changes(
            manual_df,
            args.output_dir / "manual_top_relative_change.png",
            "Manual top changes",
            args.plot_top_n,
            max_abs,
        )

    if auto_df is not None and manual_df is not None:
        _plot_top_overlap(auto_df, manual_df, args.output_dir / "top_change_overlap.png", args.plot_top_n)
        _plot_overlap_all(auto_df, manual_df, args.output_dir / "overall_overlap.png")
        _plot_relative_change_scatter(
            auto_df,
            manual_df,
            args.output_dir / "relative_change_scatter.png",
        )


if __name__ == "__main__":
    main()
