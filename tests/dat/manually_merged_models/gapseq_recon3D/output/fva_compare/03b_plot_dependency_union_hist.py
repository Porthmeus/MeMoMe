#!/usr/bin/env python3
"""
Extra diagnostic plot: dependency (relative_change) histogram for *both models separately*,
without restricting to the intersection of canonicalized IDs.

This is useful to visualize the full per-model dependency distributions, since the main
plots in 03_plot_comparisons.py are based on the intersection-only table
auto_manual_host_fva_comparison.csv.

Inputs (expected)
-----------------
From this folder:
  - auto_host_fva_baseline.csv
  - auto_host_fva_closed.csv
  - manual_host_fva_baseline.csv
  - manual_host_fva_closed.csv
Optional:
  - excluded_loop_reactions.txt (IDs to exclude, generated in 02_postprocess_fva_analysis.py)

Output
------
  results/dependency_hist_union.png
  results/dependency_hist_union.pdf
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
RESULTS_DIR = HERE / "results"


def canonicalize_rxn_id(rxn_id: str) -> str:
    base = rxn_id

    for pref in ("model1_", "model2_", "H_", "M_"):
        if base.startswith(pref):
            base = base[len(pref) :]
            break

    for pref in ("sink_model1_", "sink_model2_"):
        if base.startswith(pref):
            base = "sink_" + base[len(pref) :]
            break

    for pref in ("DM_model1_", "DM_model2_"):
        if base.startswith(pref):
            base = "DM_" + base[len(pref) :]
            break

    return base


def _load_fva(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    if "range" not in df.columns:
        df["range"] = df["maximum"] - df["minimum"]
    return df


def _maybe_filter_loops(df: pd.DataFrame, excluded_ids: set[str]) -> pd.DataFrame:
    if not excluded_ids:
        return df
    return df.drop(index=excluded_ids, errors="ignore")


def _compute_relative_change(baseline: pd.DataFrame, closed: pd.DataFrame) -> pd.Series:
    # Align to common index for this model (no cross-model intersection here)
    idx = baseline.index.intersection(closed.index)
    baseline = baseline.loc[idx]
    closed = closed.loc[idx]

    baseline_range = baseline["range"].astype(float)
    closed_range = closed["range"].astype(float)

    keep = baseline_range != 0
    baseline_range = baseline_range[keep]
    closed_range = closed_range[keep]

    rel = (baseline_range - closed_range) / baseline_range
    rel = rel.replace([np.inf, -np.inf], np.nan).dropna()
    rel = rel.round(6)

    return rel


def _dedupe_by_canonical_id(rel: pd.Series) -> pd.Series:
    # If multiple reaction IDs canonicalize to the same ID, take the mean (rare; avoids double counting).
    canon = [canonicalize_rxn_id(r) for r in rel.index.astype(str)]
    # Drop connector reactions (handled separately in exchange analysis).
    keep_mask = [not c.startswith(("IEX_", "TR_")) for c in canon]
    canon = [c for c, k in zip(canon, keep_mask, strict=False) if k]
    rel = rel.iloc[[i for i, k in enumerate(keep_mask) if k]]
    tmp = pd.DataFrame({"canonical_id": canon, "relative_change": rel.values})
    out = tmp.groupby("canonical_id", sort=False)["relative_change"].mean()
    return out


def _plot_hist(auto_rel: pd.Series, manual_rel: pd.Series) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        print("matplotlib not available; skipping plot generation.")
        return

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    combined = pd.concat([auto_rel, manual_rel], axis=0)
    if combined.empty:
        raise ValueError("No dependency values available for plotting.")

    xmin = float(np.nanmin(combined.values))
    xmax = float(np.nanmax(combined.values))
    # Keep a stable, interpretable axis for dependency (expected in [0, 1])
    xmin = min(-0.05, xmin)
    xmax = max(1.05, xmax)

    bins = np.linspace(xmin, xmax, 60)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.hist(auto_rel.values, bins=bins, alpha=0.55, label=f"MeMoMe merge (n={len(auto_rel)})")
    ax.hist(manual_rel.values, bins=bins, alpha=0.55, label=f"Manual merge (n={len(manual_rel)})")
    ax.axvline(0.8, color="red", linewidth=2, linestyle="--", label="Dependency = 0.8")

    ax.set_xlabel("Dependency (relative change in FVA range)")
    ax.set_ylabel("Count")
    ax.set_title("Dependency distributions (per-model; not restricted to intersection)")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "dependency_hist_union.png", dpi=200)
    fig.savefig(RESULTS_DIR / "dependency_hist_union.pdf")
    plt.close(fig)


def main() -> None:
    excluded_loop_path = HERE / "excluded_loop_reactions.txt"
    excluded_ids: set[str] = set()
    if excluded_loop_path.exists():
        excluded_ids = {line.strip() for line in excluded_loop_path.read_text().splitlines() if line.strip()}

    auto_baseline = _load_fva(HERE / "auto_host_fva_baseline.csv")
    auto_closed = _load_fva(HERE / "auto_host_fva_closed.csv")
    manual_baseline = _load_fva(HERE / "manual_host_fva_baseline.csv")
    manual_closed = _load_fva(HERE / "manual_host_fva_closed.csv")

    auto_baseline = _maybe_filter_loops(auto_baseline, excluded_ids)
    auto_closed = _maybe_filter_loops(auto_closed, excluded_ids)
    manual_baseline = _maybe_filter_loops(manual_baseline, excluded_ids)
    manual_closed = _maybe_filter_loops(manual_closed, excluded_ids)

    auto_rel = _dedupe_by_canonical_id(_compute_relative_change(auto_baseline, auto_closed))
    manual_rel = _dedupe_by_canonical_id(_compute_relative_change(manual_baseline, manual_closed))

    _plot_hist(auto_rel=auto_rel, manual_rel=manual_rel)

    print(
        "Wrote:",
        RESULTS_DIR / "dependency_hist_union.png",
        "and",
        RESULTS_DIR / "dependency_hist_union.pdf",
    )


if __name__ == "__main__":
    main()
