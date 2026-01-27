#!/usr/bin/env python3
"""Analyze microbiome-dependency of host reactions and compare auto vs manual models."""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import pandas as pd


LOOP_FLUX_CUTOFF = 999.0
DEPENDENCY_THRESHOLD = 0.5
FULL_DEPENDENCY_THRESHOLD = -0.999999


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


def canonical_reaction_id(rxn_id: str) -> str:
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


def _load_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df = df.round(6)
    if "minimum" in df.columns and "maximum" in df.columns:
        df["range"] = df["maximum"] - df["minimum"]
    return df


def _canonical_rel_change(df: pd.DataFrame) -> pd.Series:
    rel = df["relative_change"].dropna()
    grouped: dict[str, list[float]] = {}
    for rid, value in rel.items():
        grouped.setdefault(canonical_reaction_id(rid), []).append(float(value))
    return pd.Series({rid: sum(vals) / len(vals) for rid, vals in grouped.items()})


def _loop_reactions(baseline_df: pd.DataFrame, cutoff: float) -> set[str]:
    looped = baseline_df[
        (baseline_df["minimum"].abs() >= cutoff) | (baseline_df["maximum"].abs() >= cutoff)
    ]
    return {canonical_reaction_id(rid) for rid in looped.index}


def _write_set(path: Path, values: set[str]) -> None:
    pd.Series(sorted(values), name="reaction_id").to_csv(path, index=False)


def _write_series(path: Path, series: pd.Series, column: str) -> None:
    df = series.sort_values().to_frame(name=column)
    df.to_csv(path)


def _plot_delta_scatter(
    scatter_df: pd.DataFrame,
    top_delta: pd.Series,
    names: pd.Series,
    output_path: Path,
    top_n: int = 10,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping scatter plot.")
        return

    highlight = set(top_delta.index)
    xs = scatter_df["auto_rel_change"]
    ys = scatter_df["manual_rel_change"]
    colors = ["#D62728" if idx in highlight else "#4C72B0" for idx in scatter_df.index]

    max_abs = max(xs.abs().max(), ys.abs().max(), 1.0)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(xs, ys, c=colors, s=12, alpha=0.7, edgecolors="none")
    ax.plot([-max_abs, max_abs], [-max_abs, max_abs], linestyle="--", color="#999999", linewidth=1)
    ax.set_xlim(-max_abs, max_abs)
    ax.set_ylim(-max_abs, max_abs)
    ax.set_xlabel("Auto relative change")
    ax.set_ylabel("Manual relative change")
    ax.set_title("Relative change per reaction (filtered)")

    for rid in top_delta.head(top_n).index:
        if rid not in scatter_df.index:
            continue
        x = scatter_df.loc[rid, "auto_rel_change"]
        y = scatter_df.loc[rid, "manual_rel_change"]
        name = names.get(rid, "")
        label = f"{rid}" if not name else f"{rid}: {name}"
        ax.annotate(label, (x, y), textcoords="offset points", xytext=(6, 4), fontsize=8)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare microbiome dependency between auto and manual models.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare"),
    )
    parser.add_argument(
        "--loop-cutoff",
        type=float,
        default=LOOP_FLUX_CUTOFF,
        help="Absolute flux threshold to mark loop reactions in baseline FVA.",
    )
    parser.add_argument(
        "--dependency-threshold",
        type=float,
        default=DEPENDENCY_THRESHOLD,
        help="Dependency threshold for relative_change (absolute value).",
    )
    parser.add_argument(
        "--exclude-loops",
        choices=["none", "union", "intersection"],
        default="union",
        help="Which loop set to exclude from dependency comparison.",
    )
    parser.add_argument(
        "--top-delta",
        type=int,
        default=50,
        help="Number of reactions with largest abs(delta) to report.",
    )
    parser.add_argument(
        "--label-top-n",
        type=int,
        default=10,
        help="Number of top delta reactions to label on the scatter plot.",
    )
    args = parser.parse_args()

    auto_baseline = _load_csv(args.output_dir / "auto_host_fva_baseline.csv")
    manual_baseline = _load_csv(args.output_dir / "manual_host_fva_baseline.csv")
    auto_comp = _load_csv(args.output_dir / "auto_host_fva_comparison.csv")
    manual_comp = _load_csv(args.output_dir / "manual_host_fva_comparison.csv")

    auto_loop = _loop_reactions(auto_baseline, args.loop_cutoff)
    manual_loop = _loop_reactions(manual_baseline, args.loop_cutoff)

    loop_dir = args.output_dir / "microbiome_dependency"
    loop_dir.mkdir(parents=True, exist_ok=True)

    only_auto = auto_loop - manual_loop
    only_manual = manual_loop - auto_loop
    overlap = auto_loop & manual_loop
    union_minus_intersection = (auto_loop | manual_loop) - overlap
    _write_set(loop_dir / "loop_only_auto.csv", only_auto)
    _write_set(loop_dir / "loop_only_manual.csv", only_manual)
    _write_set(loop_dir / "loop_overlap.csv", overlap)

    loop_set = set()
    if args.exclude_loops == "union":
        loop_set = auto_loop | manual_loop
    elif args.exclude_loops == "intersection":
        loop_set = auto_loop & manual_loop

    auto_rel = _canonical_rel_change(auto_comp)
    manual_rel = _canonical_rel_change(manual_comp)
    overlap = sorted(set(auto_rel.index) & set(manual_rel.index))

    filtered = [rid for rid in overlap if rid not in loop_set]
    auto_filtered = auto_rel.loc[filtered]
    manual_filtered = manual_rel.loc[filtered]

    scatter_df = pd.DataFrame(
        {
            "auto_rel_change": auto_filtered,
            "manual_rel_change": manual_filtered,
        }
    )
    scatter_df.to_csv(loop_dir / "relative_change_scatter_data.csv")

    dependency_cutoff = -abs(args.dependency_threshold)
    auto_dependent = auto_filtered[auto_filtered <= dependency_cutoff]
    manual_dependent = manual_filtered[manual_filtered <= dependency_cutoff]
    auto_full = auto_filtered[auto_filtered <= FULL_DEPENDENCY_THRESHOLD]
    manual_full = manual_filtered[manual_filtered <= FULL_DEPENDENCY_THRESHOLD]

    _write_series(loop_dir / "auto_high_dependency.csv", auto_dependent, "auto_rel_change")
    _write_series(loop_dir / "manual_high_dependency.csv", manual_dependent, "manual_rel_change")
    _write_series(loop_dir / "auto_full_dependency.csv", auto_full, "auto_rel_change")
    _write_series(loop_dir / "manual_full_dependency.csv", manual_full, "manual_rel_change")

    delta = (manual_filtered - auto_filtered).abs().sort_values(ascending=False)
    top_delta = delta.head(args.top_delta)
    top_delta.to_frame(name="abs_delta_rel_change").to_csv(loop_dir / "top_delta_rel_change.csv")

    name_series = pd.Series(dtype=str)
    if "name" in auto_comp.columns:
        name_series = auto_comp["name"].dropna()
    if "name" in manual_comp.columns:
        manual_names = manual_comp["name"].dropna()
        name_series = name_series.combine_first(manual_names)
    name_series = name_series.groupby(name_series.index).first()
    canonical_names = {
        canonical_reaction_id(rid): name for rid, name in name_series.items() if name
    }
    names = pd.Series(canonical_names)
    _plot_delta_scatter(
        scatter_df,
        top_delta,
        names,
        loop_dir / "relative_change_scatter_highlight.png",
        top_n=args.label_top_n,
    )

    summary_lines = [
        f"loop_cutoff: {args.loop_cutoff}",
        f"loop_auto: {len(auto_loop)}",
        f"loop_manual: {len(manual_loop)}",
        f"loop_overlap: {len(overlap)}",
        f"loop_union_minus_intersection: {len(union_minus_intersection)}",
        f"exclude_loops: {args.exclude_loops}",
        f"remaining_overlap: {len(filtered)}",
        f"auto_high_dependency: {len(auto_dependent)}",
        f"manual_high_dependency: {len(manual_dependent)}",
        f"auto_full_dependency: {len(auto_full)}",
        f"manual_full_dependency: {len(manual_full)}",
    ]
    (loop_dir / "summary.txt").write_text("\n".join(summary_lines))

    print("Remaining overlap:", len(filtered))
    print("Loop union minus intersection:", len(union_minus_intersection))


if __name__ == "__main__":
    main()
