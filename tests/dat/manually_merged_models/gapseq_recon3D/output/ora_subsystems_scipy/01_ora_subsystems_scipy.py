#!/usr/bin/env python3
"""
Reimplementation of the ORA analysis (04_ora_subsystems.py) using SciPy.

This script mirrors the original workflow:
1) Define high-dependency reaction sets from dependency values.
2) Load subsystem annotations from the merged SBML (auto model) and map them to canonical IDs.
3) Run ORA with an upper-tail hypergeometric test and BH-FDR correction.

Key difference: the hypergeometric p-values and BH adjustment are computed with SciPy:
  - scipy.stats.hypergeom.sf
  - scipy.stats.false_discovery_control(method="bh")

Outputs are written to:
  tests/dat/manually_merged_models/gapseq_recon3D/output/ora_subsystems_scipy/
and compared against the reference outputs in:
  tests/dat/manually_merged_models/gapseq_recon3D/output/ora_subsystems/results_*
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, hypergeom


ROOT = Path("tests/dat/manually_merged_models/gapseq_recon3D/output")
FVA_DIR = ROOT / "fva_compare"

AUTO_MODEL_PATH = ROOT / "automatically_merged_metamodel.xml"
POSTPROCESS_COMPARISON_PATH = FVA_DIR / "auto_manual_host_fva_comparison.csv"
AUTO_BASELINE_PATH = FVA_DIR / "auto_host_fva_baseline.csv"
AUTO_CLOSED_PATH = FVA_DIR / "auto_host_fva_closed.csv"
MANUAL_BASELINE_PATH = FVA_DIR / "manual_host_fva_baseline.csv"
MANUAL_CLOSED_PATH = FVA_DIR / "manual_host_fva_closed.csv"

REFERENCE_EXCLUDE = ROOT / "ora_subsystems" / "results_exclude_loops"
REFERENCE_INCLUDE = ROOT / "ora_subsystems" / "results_include_loops"


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
        out.setdefault(cid, sub)
    return out


def load_background_and_subsystems(
    *,
    canonical_ids: set[str],
    drop_subsystems: list[str],
    fallback_background_csv: Path | None,
) -> tuple[set[str], dict[str, str]]:
    """
    Return (background_ids, subsystem_by_id).

    Preferred: load subsystems from the auto merged SBML via cobrapy.
    Fallback (when cobra is unavailable): load a precomputed background mapping from CSV.
    """
    dropped = set(drop_subsystems)
    try:
        subsystem_by_id = load_subsystem_map_from_model(AUTO_MODEL_PATH, canonical_ids)
        background_ids = {rid for rid, sub in subsystem_by_id.items() if sub not in dropped}
        return background_ids, subsystem_by_id
    except ModuleNotFoundError:
        if fallback_background_csv is None:
            raise
        bg = pd.read_csv(fallback_background_csv)
        if "reaction_id" not in bg.columns or "subsystem" not in bg.columns:
            raise ValueError(f"{fallback_background_csv} missing required columns")
        bg["reaction_id"] = bg["reaction_id"].astype(str)
        bg["subsystem"] = bg["subsystem"].astype(str)
        subsystem_by_id = dict(zip(bg["reaction_id"], bg["subsystem"], strict=True))
        background_ids = {rid for rid, sub in subsystem_by_id.items() if sub not in dropped}
        print(
            f"[warn] cobra not available; loaded background+subsystems from {fallback_background_csv}"
        )
        return background_ids, subsystem_by_id


def run_ora_scipy(
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

        # X ~ Hypergeom(N, K, n). Upper-tail p-value P[X >= k] = sf(k-1).
        p = float(hypergeom.sf(k - 1, N, K, n))
        enrich = ((k / n) / (K / N)) if (n > 0 and K > 0) else np.nan
        rows.append(
            {
                "subsystem": sub,
                "k_high_in_subsystem": int(k),
                "n_high_total": int(n),
                "K_background_in_subsystem": int(K),
                "N_background_total": int(N),
                "enrichment": float(enrich),
                "p_value": p,
            }
        )

    res = pd.DataFrame(rows)
    if res.empty:
        return res

    res["fdr_bh"] = false_discovery_control(res["p_value"].to_numpy(dtype=float), method="bh")
    res = res.sort_values(["fdr_bh", "p_value", "enrichment"], ascending=[True, True, False])
    return res


def _write_background(output_dir: Path, background_ids: set[str], subsystem_by_id: dict[str, str]) -> None:
    bg_df = pd.DataFrame(
        {
            "reaction_id": sorted(background_ids),
            "subsystem": [subsystem_by_id[r] for r in sorted(background_ids)],
        }
    )
    bg_df.to_csv(output_dir / "ora_background.csv", index=False)


def _write_high_dependency(
    *,
    output_dir: Path,
    model_label: str,
    thr: float,
    dep_bg: pd.Series,
    subsystem_by_id: dict[str, str],
) -> set[str]:
    high_ids = set(dep_bg.index[dep_bg >= thr].tolist())
    high_df = pd.DataFrame(
        {
            "reaction_id": sorted(high_ids),
            "dependency": [float(dep_bg.loc[r]) for r in sorted(high_ids)],
            "subsystem": [subsystem_by_id.get(r, "") for r in sorted(high_ids)],
        }
    )
    high_df.to_csv(output_dir / f"{model_label}_high_dependency_ge_{thr:.1f}.csv", index=False)
    return high_ids


def _compare_csv_tables(*, new_path: Path, ref_path: Path, key: str) -> dict[str, object]:
    if not ref_path.exists():
        return {"file": new_path.name, "status": "missing_reference"}
    if not new_path.exists():
        return {"file": new_path.name, "status": "missing_new"}

    new = pd.read_csv(new_path)
    ref = pd.read_csv(ref_path)

    out: dict[str, object] = {"file": new_path.name, "status": "ok"}

    if key not in new.columns or key not in ref.columns:
        out["status"] = "missing_key"
        return out

    new_key = set(new[key].astype(str).tolist())
    ref_key = set(ref[key].astype(str).tolist())
    out["n_rows_new"] = int(len(new))
    out["n_rows_ref"] = int(len(ref))
    out["n_key_new"] = int(len(new_key))
    out["n_key_ref"] = int(len(ref_key))
    out["missing_in_new"] = int(len(ref_key - new_key))
    out["missing_in_ref"] = int(len(new_key - ref_key))

    # Join on key for numeric comparisons where possible
    common = sorted(new_key & ref_key)
    if not common:
        return out

    new = new.copy()
    ref = ref.copy()
    new[key] = new[key].astype(str)
    ref[key] = ref[key].astype(str)
    new = new.set_index(key).loc[common]
    ref = ref.set_index(key).loc[common]

    for col in set(new.columns) & set(ref.columns):
        if col == key:
            continue
        # Only compare numeric columns
        try:
            a = pd.to_numeric(new[col], errors="coerce").astype(float)
            b = pd.to_numeric(ref[col], errors="coerce").astype(float)
        except Exception:
            continue
        if a.notna().any() and b.notna().any():
            diff = (a - b).abs()
            out[f"max_abs_diff__{col}"] = float(diff.max())

    return out


def compare_to_reference(*, output_dir: Path, reference_dir: Path) -> pd.DataFrame:
    files = [
        "ora_background.csv",
        "auto_high_dependency_ge_0.4.csv",
        "auto_high_dependency_ge_0.8.csv",
        "manual_high_dependency_ge_0.4.csv",
        "manual_high_dependency_ge_0.8.csv",
        "auto_ora_ge_0.4.csv",
        "auto_ora_ge_0.8.csv",
        "manual_ora_ge_0.4.csv",
        "manual_ora_ge_0.8.csv",
    ]
    rows = []
    for fname in files:
        key = "reaction_id" if "dependency" in fname or "background" in fname else "subsystem"
        rows.append(
            _compare_csv_tables(
                new_path=output_dir / fname,
                ref_path=reference_dir / fname,
                key=key,
            )
        )
    return pd.DataFrame(rows)


def run_once(*, loop_mode: str, output_dir: Path, thresholds: list[float], drop_subsystems: list[str]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    if loop_mode == "exclude":
        comp = load_postprocess_comparison(POSTPROCESS_COMPARISON_PATH)
        canonical_ids = set(comp.index.tolist())
        dep_auto = comp["relative_change_auto"].astype(float)
        dep_manual = comp["relative_change_manual"].astype(float)
    elif loop_mode == "include":
        dep_auto = compute_dependency_from_raw(baseline_path=AUTO_BASELINE_PATH, closed_path=AUTO_CLOSED_PATH)
        dep_manual = compute_dependency_from_raw(
            baseline_path=MANUAL_BASELINE_PATH, closed_path=MANUAL_CLOSED_PATH
        )
        common = dep_auto.index.intersection(dep_manual.index)
        dep_auto = dep_auto.loc[common]
        dep_manual = dep_manual.loc[common]
        canonical_ids = set(common.tolist())
    else:
        raise ValueError(f"Unknown loop_mode: {loop_mode!r}")

    fallback_ref = REFERENCE_EXCLUDE if loop_mode == "exclude" else REFERENCE_INCLUDE
    background_ids, subsystem_by_id = load_background_and_subsystems(
        canonical_ids=canonical_ids,
        drop_subsystems=drop_subsystems,
        fallback_background_csv=fallback_ref / "ora_background.csv",
    )
    _write_background(output_dir, background_ids, subsystem_by_id)

    for model_label, dep in (("auto", dep_auto), ("manual", dep_manual)):
        dep_bg = dep.loc[dep.index.intersection(pd.Index(sorted(background_ids)))]
        for thr in thresholds:
            high_ids = _write_high_dependency(
                output_dir=output_dir,
                model_label=model_label,
                thr=thr,
                dep_bg=dep_bg,
                subsystem_by_id=subsystem_by_id,
            )

            ora = run_ora_scipy(background_ids=background_ids, high_ids=high_ids, subsystem_by_id=subsystem_by_id)
            ora.to_csv(output_dir / f"{model_label}_ora_ge_{thr:.1f}.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-root",
        type=Path,
        default=ROOT / "ora_subsystems_scipy",
        help="Root output directory for SciPy-based ORA reimplementation.",
    )
    parser.add_argument(
        "--loop-modes",
        choices=["exclude", "include"],
        nargs="+",
        default=["exclude", "include"],
        help="Which loop modes to run.",
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
    parser.add_argument(
        "--compare",
        action="store_true",
        help="Compare SciPy outputs against the reference outputs produced by 04_ora_subsystems.py.",
    )
    args = parser.parse_args()

    for mode in args.loop_modes:
        out_dir = args.output_root / f"results_{mode}_loops"
        run_once(
            loop_mode=mode,
            output_dir=out_dir,
            thresholds=list(args.thresholds),
            drop_subsystems=list(args.drop_subsystems),
        )

        if args.compare:
            ref_dir = REFERENCE_EXCLUDE if mode == "exclude" else REFERENCE_INCLUDE
            report = compare_to_reference(output_dir=out_dir, reference_dir=ref_dir)
            report.to_csv(out_dir / "comparison_to_reference.csv", index=False)
            print(f"[{mode}] wrote comparison report to {out_dir / 'comparison_to_reference.csv'}")


if __name__ == "__main__":
    main()
