#!/usr/bin/env python3
"""Postprocess FVA CSVs: round, compare, and summarize."""

from __future__ import annotations

import argparse
from pathlib import Path

import cobra
import pandas as pd

FLUX_TOL = 1e-6


def _load_fva(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df = df.round(6)
    if "minimum" in df.columns and "maximum" in df.columns:
        df["range"] = df["maximum"] - df["minimum"]
    return df


def _load_names(model_path: Path) -> pd.Series:
    model = cobra.io.read_sbml_model(str(model_path))
    return pd.Series({rxn.id: (rxn.name or "") for rxn in model.reactions})


def _join_ranges(
    baseline: pd.DataFrame, closed: pd.DataFrame | None, names: pd.Series
) -> pd.DataFrame:
    baseline = baseline.rename(
        columns={"minimum": "min_before", "maximum": "max_before", "range": "range_before"}
    )
    if closed is None:
        closed_block = pd.DataFrame(
            {
                "min_after": pd.NA,
                "max_after": pd.NA,
                "range_after": pd.NA,
            },
            index=baseline.index,
        )
        merged = baseline.join(closed_block, how="left")
    else:
        closed = closed.rename(
            columns={"minimum": "min_after", "maximum": "max_after", "range": "range_after"}
        )
        merged = baseline.join(closed[["min_after", "max_after", "range_after"]], how="left")
    merged.insert(0, "name", names)

    def pct_change(row: pd.Series) -> float | None:
        before = row["range_before"]
        after = row["range_after"]
        if before == 0 and after == 0:
            return None
        if before == 0:
            return float("inf")
        return (after - before) / before * 100.0

    merged["pct_change"] = merged.apply(pct_change, axis=1)
    merged["range_delta"] = merged["range_after"] - merged["range_before"]

    def rel_change(row: pd.Series) -> float | None:
        before = row["range_before"]
        after = row["range_after"]
        denom = max(before, after)
        if denom == 0:
            return None
        return (after - before) / denom

    merged["relative_change"] = merged.apply(rel_change, axis=1)
    return merged.round(6)


def _summarize(df: pd.DataFrame) -> dict[str, float]:
    zero_before = (df["range_before"] == 0).sum()
    zero_after = (df["range_after"] == 0).sum()
    changed = (df["range_before"] != df["range_after"]).sum()
    increased = (df["range_delta"] > FLUX_TOL).sum()
    finite = df["pct_change"].replace([float("inf"), float("-inf")], pd.NA).dropna()
    rel_finite = df["relative_change"].replace([float("inf"), float("-inf")], pd.NA).dropna()
    mean_pct = float(finite.mean()) if not finite.empty else 0.0
    median_pct = float(finite.median()) if not finite.empty else 0.0
    mean_rel = float(rel_finite.mean()) if not rel_finite.empty else 0.0
    median_rel = float(rel_finite.median()) if not rel_finite.empty else 0.0
    return {
        "host_reactions": len(df),
        "zero_range_before": int(zero_before),
        "zero_range_after": int(zero_after),
        "changed_ranges": int(changed),
        "increased_ranges": int(increased),
        "mean_pct_change": mean_pct,
        "median_pct_change": median_pct,
        "mean_rel_change": mean_rel,
        "median_rel_change": median_rel,
    }


def _forced_nonzero_from_fva(fva: pd.DataFrame, output_path: Path) -> int:
    forced = fva[(fva["minimum"] > FLUX_TOL) | (fva["maximum"] < -FLUX_TOL)].copy()
    forced["range"] = forced["maximum"] - forced["minimum"]
    forced.to_csv(output_path, index=True)
    return len(forced)


def _load_closed(path: Path) -> tuple[pd.DataFrame | None, bool]:
    if not path.exists():
        return None, False
    df = pd.read_csv(path)
    if "status" in df.columns:
        status = str(df["status"].iloc[0]).lower()
        if status == "infeasible":
            return None, False
    df = pd.read_csv(path, index_col=0)
    df = df.round(6)
    if "minimum" in df.columns and "maximum" in df.columns:
        df["range"] = df["maximum"] - df["minimum"]
    return df, True


def _find_host_biomass_id(index: pd.Index) -> str | None:
    candidates = [rid for rid in index if "biomass_reaction" in str(rid)]
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]
    filtered = [rid for rid in candidates if "Bacterial" not in str(rid)]
    if len(filtered) == 1:
        return filtered[0]
    return candidates[0]


def _postprocess_one(
    label: str,
    output_dir: Path,
    model_path: Path,
) -> dict[str, float | int | bool]:
    baseline_path = output_dir / f"{label}_host_fva_baseline.csv"
    closed_path = output_dir / f"{label}_host_fva_closed.csv"
    exchange_fva_path = output_dir / f"{label}_bacterial_exchange_fva.csv"

    baseline = _load_fva(baseline_path)
    closed, closed_feasible = _load_closed(closed_path)
    names = _load_names(model_path)

    merged = _join_ranges(baseline, closed, names)
    merged.to_csv(output_dir / f"{label}_host_fva_comparison.csv")

    summary = _summarize(merged)

    biomass_id = _find_host_biomass_id(baseline.index)
    max_host_flux_open = (
        float(baseline.loc[biomass_id, "maximum"]) if biomass_id in baseline.index else float("nan")
    )
    if closed is not None and biomass_id in closed.index:
        max_host_flux_closed = float(closed.loc[biomass_id, "maximum"])
    else:
        max_host_flux_closed = float("nan")

    exchange_fva = _load_fva(exchange_fva_path) if exchange_fva_path.exists() else None
    if exchange_fva is not None and biomass_id in exchange_fva.index:
        host_only_objective = float(exchange_fva.loc[biomass_id, "maximum"])
        exchange_fva = exchange_fva.drop(index=[biomass_id])
    else:
        host_only_objective = float("nan")

    if exchange_fva is not None:
        forced_count = _forced_nonzero_from_fva(
            exchange_fva, output_dir / f"{label}_bacterial_exchange_forced_nonzero.csv"
        )
        exchange_nonzero = exchange_fva[
            (exchange_fva["minimum"].abs() > FLUX_TOL)
            | (exchange_fva["maximum"].abs() > FLUX_TOL)
        ]
        exchange_nonzero.to_csv(
            output_dir / f"{label}_bacterial_exchange_nonzero.csv", index=True
        )
        nonzero_count = len(exchange_nonzero)
        exchange_count = len(exchange_fva)

        def is_forced_zero_id(rid: str) -> bool:
            if "cpd00007" not in rid:
                return False
            return rid.startswith("TR_model2_") or "IEX_" in rid

        forced_zero_count = exchange_fva[
            exchange_fva.apply(
                lambda row: is_forced_zero_id(str(row.name))
                and row["minimum"] == 0
                and row["maximum"] == 0,
                axis=1,
            )
        ].shape[0]
    else:
        forced_count = 0
        nonzero_count = 0
        exchange_count = 0
        forced_zero_count = 0

    summary["bacterial_exchange_forced_nonzero"] = forced_count
    summary["bacterial_exchange_nonzero_open"] = nonzero_count
    summary["bacterial_exchange_closed"] = exchange_count
    summary["closed_feasible"] = bool(closed_feasible)
    summary["max_host_flux_open"] = max_host_flux_open
    summary["max_host_flux_closed"] = max_host_flux_closed
    summary["host_only_objective"] = host_only_objective
    summary["forced_zero_exchange_count"] = forced_zero_count
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Postprocess FVA CSVs for auto/manual models.")
    parser.add_argument(
        "--auto-model",
        type=Path,
        default=Path(
            "tests/dat/manually_merged_models/gapseq_recon3D/output/automatically_merged_metamodel.xml"
        ),
    )
    parser.add_argument(
        "--manual-model",
        type=Path,
        default=Path(
            "tests/dat/manually_merged_models/gapseq_recon3D/output/manual_diet/"
            "merged_model_2025_prefixed_normalized_diet.xml"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare"),
    )
    args = parser.parse_args()

    summaries = {
        "auto": _postprocess_one("auto", args.output_dir, args.auto_model),
        "manual": _postprocess_one("manual", args.output_dir, args.manual_model),
    }

    summary_path = args.output_dir / "summary.txt"
    lines = []
    for label, summary in summaries.items():
        lines.append(f"[{label}]")
        for key, value in summary.items():
            lines.append(f"{key}: {value}")
        lines.append("")
    summary_path.write_text("\n".join(lines))

    for label, summary in summaries.items():
        print(label, summary)


if __name__ == "__main__":
    main()
