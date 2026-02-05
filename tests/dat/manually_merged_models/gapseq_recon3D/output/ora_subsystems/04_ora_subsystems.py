#!/usr/bin/env python3
"""
Over-Representation Analysis (ORA) of host subsystems among "high dependency" reactions.

Dependency is computed from host FVA range widths:
    dependency = (range_open - range_closed) / range_open

This script uses the *postprocessed* comparison table produced by
02_postprocess_fva_analysis.py:
    fva_compare/auto_manual_host_fva_comparison.csv

That table already:
- canonicalizes IDs across auto/manual
- intersects the two models
- typically excludes loops/blocked reactions (depending on the postprocess settings)

Loop handling:
- `--loop-mode exclude` (default): use the postprocessed comparison table (typically loop-free).
- `--loop-mode include`: recompute dependency directly from the raw baseline/closed FVA CSVs,
  keeping loop reactions, but still dropping blocked reactions (range==0 in both baseline
  and closed) to match the postprocess behaviour.

We then:
1) Build two high-dependency sets per model using thresholds (default: 0.4 and 0.8).
2) Load subsystem annotations from the merged SBML (auto) and map them to canonical IDs.
3) Run ORA (hypergeometric upper-tail) and apply BH-FDR.

Background is restricted to reactions with a non-empty subsystem annotation, excluding
dropped subsystems (e.g. Transport/exchange).
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("tests/dat/manually_merged_models/gapseq_recon3D/output")
FVA_DIR = ROOT / "fva_compare"

AUTO_MODEL_PATH = ROOT / "automatically_merged_metamodel.xml"
MANUAL_MODEL_PATH = ROOT / "merged_model_2025_prefixed_normalized_diet_fixIEX.xml"
POSTPROCESS_COMPARISON_PATH = FVA_DIR / "auto_manual_host_fva_comparison.csv"
AUTO_BASELINE_PATH = FVA_DIR / "auto_host_fva_baseline.csv"
AUTO_CLOSED_PATH = FVA_DIR / "auto_host_fva_closed.csv"
MANUAL_BASELINE_PATH = FVA_DIR / "manual_host_fva_baseline.csv"
MANUAL_CLOSED_PATH = FVA_DIR / "manual_host_fva_closed.csv"


def canonicalize_rxn_id(rxn_id: str) -> str:
    # Canonicalize the IDs as they appear in *FVA outputs* and *cobra model*:
    # - auto host: model1_*
    # - manual host: H_*
    # - occasional prefixes inside IDs: DM_model1_*
    base = str(rxn_id)

    # Leading model prefixes
    for pref in ("model1_", "model2_", "H_", "M_"):
        if base.startswith(pref):
            base = base[len(pref) :]
            break

    # Token-before-model-prefix cases (e.g. DM_model1_foo[c] -> DM_foo[c])
    for pref in ("DM_model1_", "DM_model2_"):
        if base.startswith(pref):
            base = "DM_" + base[len(pref) :]
            break

    return base


def load_postprocess_comparison(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str)
    required = {
        "relative_change_auto",
        "relative_change_manual",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {sorted(missing)}")
    return df


def load_fva(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str)
    for col in ("minimum", "maximum", "range"):
        if col not in df.columns:
            raise ValueError(f"{path} missing column '{col}'")
    return df


def compute_dependency_from_raw(*, baseline_path: Path, closed_path: Path) -> pd.Series:
    baseline = load_fva(baseline_path)
    closed = load_fva(closed_path)

    baseline.index = [canonicalize_rxn_id(i) for i in baseline.index]
    closed.index = [canonicalize_rxn_id(i) for i in closed.index]

    common = baseline.index.intersection(closed.index)
    b = baseline.loc[common, "range"].astype(float)
    c = closed.loc[common, "range"].astype(float)

    blocked = (b == 0.0) & (c == 0.0)
    b = b.loc[~blocked]
    c = c.loc[~blocked]

    dep = pd.Series(np.nan, index=b.index, dtype=float)
    nonzero = b > 0.0
    dep.loc[nonzero] = (b.loc[nonzero] - c.loc[nonzero]) / b.loc[nonzero]
    return dep


def _log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k: int, N: int, K: int, n: int) -> float:
    """
    Hypergeometric survival function P[X >= k] for X ~ Hypergeom(N, K, n).

    Implemented without scipy using log-combinations for numerical stability.
    """
    if k <= 0:
        return 1.0
    max_k = min(K, n)
    if k > max_k:
        return 0.0

    log_denom = _log_choose(N, n)
    logs = []
    for i in range(k, max_k + 1):
        log_num = _log_choose(K, i) + _log_choose(N - K, n - i)
        logs.append(log_num - log_denom)

    m = max(logs)
    return float(math.exp(m) * sum(math.exp(x - m) for x in logs))


def benjamini_hochberg(pvals: pd.Series) -> pd.Series:
    p = pvals.astype(float)
    m = int(p.notna().sum())
    if m == 0:
        return pd.Series(np.nan, index=p.index, dtype=float)

    order = p.sort_values().index
    q = pd.Series(np.nan, index=p.index, dtype=float)

    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = m - rank + 1
        val = p.loc[idx] * m / i
        prev = min(prev, val)
        q.loc[idx] = prev
    return q.clip(lower=0.0, upper=1.0)


def load_subsystem_map_from_model(model_path: Path, canonical_ids: set[str]) -> dict[str, str]:
    import cobra

    model = cobra.io.read_sbml_model(str(model_path))

    out: dict[str, str] = {}
    for rxn in model.reactions:
        cid = canonicalize_rxn_id(rxn.id)
        if cid not in canonical_ids:
            continue
        sub = rxn.subsystem
        if sub is None:
            continue
        sub = str(sub).strip()
        if not sub or sub.lower() == "nan":
            continue
        # Prefer first seen; subsystems should be consistent.
        out.setdefault(cid, sub)

    return out


def run_ora(
    *,
    background_ids: set[str],
    high_ids: set[str],
    subsystem_by_id: dict[str, str],
) -> pd.DataFrame:
    # Build subsystem -> ids for background
    subs_to_ids: dict[str, set[str]] = {}
    for rid in background_ids:
        sub = subsystem_by_id.get(rid)
        if sub is None:
            continue
        subs_to_ids.setdefault(sub, set()).add(rid)

    N = len(background_ids)
    n = len(high_ids)

    rows = []
    for sub, ids in subs_to_ids.items():
        K = len(ids)
        k = len(ids & high_ids)
        if K == 0:
            continue

        p = hypergeom_sf(k, N=N, K=K, n=n)
        enrich = ((k / n) / (K / N)) if (n > 0 and K > 0) else np.nan
        rows.append(
            {
                "subsystem": sub,
                "k_high_in_subsystem": k,
                "n_high_total": n,
                "K_background_in_subsystem": K,
                "N_background_total": N,
                "enrichment": enrich,
                "p_value": p,
            }
        )

    res = pd.DataFrame(rows)
    if res.empty:
        return res

    res["fdr_bh"] = benjamini_hochberg(res["p_value"])
    res = res.sort_values(["fdr_bh", "p_value", "enrichment"], ascending=[True, True, False])
    return res


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=POSTPROCESS_COMPARISON_PATH,
        help="Postprocess comparison CSV (from 02_postprocess_fva_analysis.py).",
    )
    parser.add_argument(
        "--loop-mode",
        choices=["exclude", "include"],
        default="exclude",
        help="Exclude loops (use postprocess table) or include loops (recompute from raw FVA).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "ora_subsystems" / "results",
        help="Output directory for high-dependency lists and ORA result tables.",
    )
    parser.add_argument(
        "--thresholds",
        type=float,
        nargs="+",
        default=[0.4, 0.8],
        help="Dependency thresholds for high-dependency sets.",
    )
    parser.add_argument(
        "--drop-subsystems",
        type=str,
        nargs="*",
        default=["Transport, extracellular", "Exchange/demand reaction"],
        help="Subsystems to drop from background and ORA.",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.loop_mode == "exclude":
        comp = load_postprocess_comparison(args.input)
        canonical_ids = set(comp.index.tolist())
        dep_auto = comp["relative_change_auto"].astype(float)
        dep_manual = comp["relative_change_manual"].astype(float)
    else:
        dep_auto = compute_dependency_from_raw(baseline_path=AUTO_BASELINE_PATH, closed_path=AUTO_CLOSED_PATH)
        dep_manual = compute_dependency_from_raw(baseline_path=MANUAL_BASELINE_PATH, closed_path=MANUAL_CLOSED_PATH)
        common = dep_auto.index.intersection(dep_manual.index)
        dep_auto = dep_auto.loc[common]
        dep_manual = dep_manual.loc[common]
        canonical_ids = set(common.tolist())

    # Subsystems from auto model (host namespace). This should be stable across both.
    subsystem_by_id = load_subsystem_map_from_model(AUTO_MODEL_PATH, canonical_ids)

    # Build background: canonical IDs with a non-dropped subsystem.
    dropped = set(args.drop_subsystems)
    background_ids = {rid for rid, sub in subsystem_by_id.items() if sub not in dropped}

    # Export background
    bg_df = pd.DataFrame(
        {
            "reaction_id": sorted(background_ids),
            "subsystem": [subsystem_by_id[r] for r in sorted(background_ids)],
        }
    )
    bg_df.to_csv(args.output_dir / "ora_background.csv", index=False)

    # ORA per model + threshold
    for model_label, dep in (("auto", dep_auto), ("manual", dep_manual)):
        # Restrict dependency vector to background ids
        dep_bg = dep.loc[dep.index.intersection(pd.Index(sorted(background_ids)))]

        for thr in args.thresholds:
            high_ids = set(dep_bg.index[dep_bg >= thr].tolist())

            high_df = pd.DataFrame(
                {
                    "reaction_id": sorted(high_ids),
                    "dependency": [float(dep_bg.loc[r]) for r in sorted(high_ids)],
                    "subsystem": [subsystem_by_id.get(r, "") for r in sorted(high_ids)],
                }
            )
            high_path = args.output_dir / f"{model_label}_high_dependency_ge_{thr:.1f}.csv"
            high_df.to_csv(high_path, index=False)

            ora = run_ora(background_ids=background_ids, high_ids=high_ids, subsystem_by_id=subsystem_by_id)
            ora_path = args.output_dir / f"{model_label}_ora_ge_{thr:.1f}.csv"
            ora.to_csv(ora_path, index=False)

            print(
                f"{model_label} thr>={thr:.1f}: background={len(background_ids)} "
                f"high={len(high_ids)} subsystems={len(ora)}"
            )


if __name__ == "__main__":
    main()
