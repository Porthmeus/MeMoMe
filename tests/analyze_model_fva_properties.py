#!/usr/bin/env python3
"""Run FVA-based host reaction sensitivity analysis for auto/manual models."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import cobra
import os
import pandas as pd
from cobra.exceptions import Infeasible
from cobra.flux_analysis import flux_variability_analysis

SPECIAL_PREFIXES = ("EX_", "DM_", "SK_", "sink_")
FLUX_TOL = 1e-6


@dataclass(frozen=True)
class ModelSpec:
    label: str
    path: Path
    host_prefix: str
    micro_prefix: str
    auto_model: Path | None = None


def _find_unique_reaction(model: cobra.Model, *, prefix: str | None, key: str) -> cobra.Reaction:
    matches = [
        rxn
        for rxn in model.reactions
        if key in rxn.id and (prefix is None or rxn.id.startswith(prefix))
    ]
    if not matches:
        matches = [
            rxn
            for rxn in model.reactions
            if rxn.name and key in rxn.name and (prefix is None or rxn.id.startswith(prefix))
        ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected 1 reaction for key '{key}' prefix '{prefix}', found {len(matches)}: "
            f"{[rxn.id for rxn in matches]}"
        )
    return matches[0]


def _identify_biomass_ids(model: cobra.Model, host_prefix: str, micro_prefix: str) -> dict[str, str]:
    host = _find_unique_reaction(model, prefix=host_prefix, key="biomass_reaction")
    gram_neg = _find_unique_reaction(
        model, prefix=micro_prefix, key="Bacterial_Gram_negative_biomass_reaction"
    )
    gram_pos = _find_unique_reaction(
        model, prefix=micro_prefix, key="Bacterial_Gram_positive_biomass_reaction"
    )
    return {"host": host.id, "gram_neg": gram_neg.id, "gram_pos": gram_pos.id}


def _is_host_reaction_auto(rxn_id: str) -> bool:
    if rxn_id.startswith("TR_"):
        return False
    if rxn_id.startswith("model1_"):
        return True
    for prefix in SPECIAL_PREFIXES:
        if rxn_id.startswith(f"{prefix}model1_"):
            return True
    return False


def _is_host_reaction_manual(rxn_id: str) -> bool:
    return rxn_id.startswith("H_")


def _host_reaction_ids(model: cobra.Model, spec: ModelSpec, max_reactions: int | None) -> list[str]:
    if spec.label == "auto":
        predicate = _is_host_reaction_auto
    else:
        predicate = _is_host_reaction_manual
    host_ids = [rxn.id for rxn in model.reactions if predicate(rxn.id)]
    if max_reactions is not None:
        return host_ids[:max_reactions]
    return host_ids


def _bacterial_exchange_ids(model: cobra.Model, spec: ModelSpec) -> list[str]:
    if spec.label == "auto":
        return [rxn.id for rxn in model.reactions if rxn.id.startswith("TR_model2_")]
    bacterial_iex = []
    for rxn in model.reactions:
        if "IEX_" not in rxn.id:
            continue
        if any(met.id.startswith("M_") for met in rxn.metabolites):
            bacterial_iex.append(rxn.id)
    return bacterial_iex


def _set_objective(model: cobra.Model, rxn_id: str) -> None:
    if rxn_id in model.reactions:
        model.objective = rxn_id


def _forced_zero_exchange_ids(model: cobra.Model, spec: ModelSpec) -> list[str]:
    if spec.label == "auto":
        return [rxn.id for rxn in model.reactions if "TR_model2_cpd00007" in rxn.id]
    forced = []
    for rxn in model.reactions:
        if "IEX_" not in rxn.id or "cpd00007" not in rxn.id:
            continue
        if any(met.id.startswith("M_") for met in rxn.metabolites):
            forced.append(rxn.id)
    return forced


def _canonical_met_id(met_id: str) -> str:
    base = met_id
    for prefix in ("H_", "M_", "model1_", "model2_"):
        if base.startswith(prefix):
            base = base[len(prefix):]
            break
    if base.endswith("_t"):
        base = base[:-2]
    for suffix in ("[e]", "[c]", "[p]", "[b]", "[f]"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    for suffix in ("_e", "_c", "_p", "_b", "_f"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    return base


def _auto_ex_t_bounds(auto_model_path: Path) -> dict[str, tuple[float, float]]:
    auto_model = cobra.io.read_sbml_model(str(auto_model_path))
    bounds = {}
    for rxn in auto_model.reactions:
        if not (rxn.id.startswith("EX_") and rxn.id.endswith("_t")):
            continue
        if not rxn.metabolites:
            continue
        met = next(iter(rxn.metabolites))
        base = _canonical_met_id(met.id)
        bounds[base] = (rxn.lower_bound, rxn.upper_bound)
    return bounds


def _apply_manual_diet_from_auto(
    model: cobra.Model, spec: ModelSpec, output_dir: Path
) -> dict[str, int]:
    if spec.label != "manual" or spec.auto_model is None:
        return {"manual_diet_applied": 0, "manual_diet_candidates": 0}
    auto_bounds = _auto_ex_t_bounds(spec.auto_model)
    candidates = 0
    applied = 0
    for rxn in model.reactions:
        if not rxn.id.startswith("H_EX_"):
            continue
        candidates += 1
        host_met = next((met for met in rxn.metabolites if met.id.startswith("H_")), None)
        if host_met is None and rxn.metabolites:
            host_met = next(iter(rxn.metabolites))
        if host_met is None:
            continue
        base = _canonical_met_id(host_met.id)
        if base not in auto_bounds:
            continue
        lb, ub = auto_bounds[base]
        rxn.lower_bound = lb
        rxn.upper_bound = ub
        applied += 1

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "manual_diet_applied.txt").write_text(
        f"applied {applied} of {candidates} H_EX bounds from {spec.auto_model}\\n"
    )
    return {"manual_diet_applied": applied, "manual_diet_candidates": candidates}


def _apply_biomass_bounds(
    model: cobra.Model,
    biomass_ids: dict[str, str],
    host_lb: float,
    gram_neg_lb: float | None,
    gram_pos_lb: float | None,
    original_lbs: dict[str, float],
) -> None:
    host_rxn = model.reactions.get_by_id(biomass_ids["host"])
    gram_neg_rxn = model.reactions.get_by_id(biomass_ids["gram_neg"])
    gram_pos_rxn = model.reactions.get_by_id(biomass_ids["gram_pos"])

    host_rxn.lower_bound = host_lb
    if gram_neg_lb is None:
        gram_neg_rxn.lower_bound = original_lbs[gram_neg_rxn.id]
    else:
        gram_neg_rxn.lower_bound = gram_neg_lb
    if gram_pos_lb is None:
        gram_pos_rxn.lower_bound = original_lbs[gram_pos_rxn.id]
    else:
        gram_pos_rxn.lower_bound = gram_pos_lb


def _close_reactions(model: cobra.Model, reaction_ids: Iterable[str]) -> None:
    for rxn_id in reaction_ids:
        rxn = model.reactions.get_by_id(rxn_id)
        rxn.lower_bound = 0.0
        rxn.upper_bound = 0.0


def _run_fva(model: cobra.Model, reaction_ids: list[str], processes: int) -> pd.DataFrame:
    fva = flux_variability_analysis(
        model, reaction_list=reaction_ids, fraction_of_optimum=0.0, processes=processes
    )
    fva["range"] = fva["maximum"] - fva["minimum"]
    return fva


def _join_ranges(
    baseline: pd.DataFrame, closed: pd.DataFrame, names: pd.Series
) -> pd.DataFrame:
    baseline = baseline.rename(columns={"minimum": "min_before", "maximum": "max_before", "range": "range_before"})
    closed = closed.rename(columns={"minimum": "min_after", "maximum": "max_after", "range": "range_after"})
    merged = baseline.join(closed[["min_after", "max_after", "range_after"]], how="left")
    merged.insert(0, "name", names)

    def pct_change(row: pd.Series) -> float | None:
        before = row["range_before"]
        after = row["range_after"]
        if before == 0:
            return 0.0 if after == 0 else float("inf")
        return (after - before) / before * 100.0

    merged["pct_change"] = merged.apply(pct_change, axis=1)
    return merged


def _summarize(df: pd.DataFrame) -> dict[str, float]:
    zero_before = (df["range_before"] == 0).sum()
    zero_after = (df["range_after"] == 0).sum()
    changed = (df["range_before"] != df["range_after"]).sum()
    finite = df["pct_change"].replace([float("inf"), float("-inf")], pd.NA).dropna()
    mean_pct = float(finite.mean()) if not finite.empty else 0.0
    median_pct = float(finite.median()) if not finite.empty else 0.0
    return {
        "host_reactions": len(df),
        "zero_range_before": int(zero_before),
        "zero_range_after": int(zero_after),
        "changed_ranges": int(changed),
        "mean_pct_change": mean_pct,
        "median_pct_change": median_pct,
    }


def _diagnose_bacterial_exchange(
    model: cobra.Model,
    bacterial_exchange_ids: list[str],
    output_dir: Path,
    label: str,
) -> tuple[float, int]:
    solution = model.optimize()
    if solution.status != "optimal":
        return float("nan"), 0
    nonzero = []
    for rxn_id in bacterial_exchange_ids:
        flux = solution.fluxes.get(rxn_id, 0.0)
        if abs(flux) > FLUX_TOL:
            rxn = model.reactions.get_by_id(rxn_id)
            nonzero.append({"reaction_id": rxn_id, "name": rxn.name or "", "flux": flux})

    if nonzero:
        df = pd.DataFrame(nonzero).sort_values("flux", key=lambda s: s.abs(), ascending=False)
    else:
        df = pd.DataFrame(columns=["reaction_id", "name", "flux"])
    df.to_csv(output_dir / f"{label}_bacterial_exchange_nonzero.csv", index=False)

    return float(solution.objective_value), len(nonzero)


def _bacterial_exchange_forced_nonzero(
    model: cobra.Model,
    bacterial_exchange_ids: list[str],
    output_dir: Path,
    label: str,
    processes: int,
) -> tuple[float, int]:
    solution = model.optimize()
    if solution.status != "optimal":
        return float("nan"), 0
    fva = _run_fva(model, bacterial_exchange_ids, processes)
    forced = fva[(fva["minimum"] > FLUX_TOL) | (fva["maximum"] < -FLUX_TOL)].copy()
    forced["range"] = forced["maximum"] - forced["minimum"]
    forced.to_csv(output_dir / f"{label}_bacterial_exchange_forced_nonzero.csv")
    return float(solution.objective_value), len(forced)


def analyze_model(
    spec: ModelSpec, output_dir: Path, processes: int, max_reactions: int | None
) -> dict[str, float]:
    model = cobra.io.read_sbml_model(str(spec.path))

    host_reaction_ids = _host_reaction_ids(model, spec, max_reactions)
    names = pd.Series({rxn.id: (rxn.name or "") for rxn in model.reactions})

    biomass_ids = _identify_biomass_ids(model, spec.host_prefix, spec.micro_prefix)
    original_lbs = {
        key: model.reactions.get_by_id(rxn_id).lower_bound
        for key, rxn_id in biomass_ids.items()
    }
    output_dir.mkdir(parents=True, exist_ok=True)

    diet_summary = _apply_manual_diet_from_auto(model, spec, output_dir)
    forced_zero_ids = _forced_zero_exchange_ids(model, spec)
    baseline_model = model.copy()
    _set_objective(baseline_model, biomass_ids["host"])
    _apply_biomass_bounds(
        baseline_model,
        biomass_ids,
        host_lb=0.0001,
        gram_neg_lb=0.001,
        gram_pos_lb=0.001,
        original_lbs=original_lbs,
    )
    _close_reactions(baseline_model, forced_zero_ids)
    baseline_solution = baseline_model.optimize()
    if baseline_solution.status != "optimal":
        raise Infeasible(f"{spec.label} baseline model infeasible")
    baseline_fva = _run_fva(baseline_model, host_reaction_ids, processes)

    host_only_model = model.copy()
    _set_objective(host_only_model, biomass_ids["host"])
    _apply_biomass_bounds(
        host_only_model,
        biomass_ids,
        host_lb=0.0001,
        gram_neg_lb=0.0,
        gram_pos_lb=0.0,
        original_lbs=original_lbs,
    )
    _close_reactions(host_only_model, forced_zero_ids)
    bacterial_exchange_ids = _bacterial_exchange_ids(host_only_model, spec)
    host_only_flux, forced_nonzero_count = _bacterial_exchange_forced_nonzero(
        host_only_model, bacterial_exchange_ids, output_dir, spec.label, processes
    )

    closed_model = model.copy()
    _set_objective(closed_model, biomass_ids["host"])
    _apply_biomass_bounds(
        closed_model,
        biomass_ids,
        host_lb=0.0001,
        gram_neg_lb=0.0,
        gram_pos_lb=0.0,
        original_lbs=original_lbs,
    )
    _close_reactions(closed_model, forced_zero_ids)
    bacterial_exchange_ids = _bacterial_exchange_ids(closed_model, spec)
    _close_reactions(closed_model, bacterial_exchange_ids)
    closed_solution = closed_model.optimize()
    if closed_solution.status != "optimal":
        baseline_fva.to_csv(output_dir / f"{spec.label}_host_fva_baseline.csv")
        max_host_open, nonzero_count = _diagnose_bacterial_exchange(
            baseline_model, bacterial_exchange_ids, output_dir, spec.label
        )
        summary = {
            "host_reactions": len(host_reaction_ids),
            "bacterial_exchange_closed": len(bacterial_exchange_ids),
            "closed_feasible": False,
            "max_host_flux_open": max_host_open,
            "max_host_flux_closed": float("nan"),
            "bacterial_exchange_nonzero_open": nonzero_count,
            "bacterial_exchange_forced_nonzero": forced_nonzero_count,
            "host_only_objective": host_only_flux,
            "forced_zero_exchange_count": len(forced_zero_ids),
            **diet_summary,
        }
        return summary

    closed_fva = _run_fva(closed_model, host_reaction_ids, processes)

    merged = _join_ranges(baseline_fva, closed_fva, names)

    baseline_fva.to_csv(output_dir / f"{spec.label}_host_fva_baseline.csv")
    closed_fva.to_csv(output_dir / f"{spec.label}_host_fva_closed.csv")
    merged.to_csv(output_dir / f"{spec.label}_host_fva_comparison.csv")

    max_host_open, nonzero_count = _diagnose_bacterial_exchange(
        baseline_model, bacterial_exchange_ids, output_dir, spec.label
    )
    max_host_closed = float(closed_solution.objective_value)
    summary = _summarize(merged)
    summary["bacterial_exchange_closed"] = len(bacterial_exchange_ids)
    summary["closed_feasible"] = True
    summary["max_host_flux_open"] = max_host_open
    summary["max_host_flux_closed"] = max_host_closed
    summary["bacterial_exchange_nonzero_open"] = nonzero_count
    summary["bacterial_exchange_forced_nonzero"] = forced_nonzero_count
    summary["host_only_objective"] = host_only_flux
    summary["forced_zero_exchange_count"] = len(forced_zero_ids)
    summary.update(diet_summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare FVA ranges for auto/manual models.")
    parser.add_argument(
        "--auto-path",
        type=Path,
        default=Path("tests/dat/manually_merged_models/gapseq_recon3D/output/automatically_merged_metamodel.xml"),
    )
    parser.add_argument(
        "--manual-path",
        type=Path,
        default=Path("tests/dat/manually_merged_models/gapseq_recon3D/merged_model_2025_prefixed_normalized.xml"),
    )
    parser.add_argument(
        "--auto-model",
        type=Path,
        default=Path("tests/dat/manually_merged_models/gapseq_recon3D/output/automatically_merged_metamodel.xml"),
        help="Auto model used to transfer EX_t diet to manual H_EX reactions.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("tests/dat/manually_merged_models/gapseq_recon3D/output/fva_compare"),
    )
    parser.add_argument(
        "--max-reactions",
        type=int,
        default=None,
        help="Optional cap on host reactions (debugging only).",
    )
    args = parser.parse_args()

    specs = [
        ModelSpec("auto", args.auto_path, "model1_", "model2_"),
        ModelSpec("manual", args.manual_path, "H_", "M_", args.auto_model),
    ]
    processes = max(1, os.cpu_count() or 1)

    summaries = {}
    for spec in specs:
        summaries[spec.label] = analyze_model(
            spec, args.output_dir, processes, args.max_reactions
        )

    summary_path = args.output_dir / "summary.txt"
    lines = []
    for label, summary in summaries.items():
        lines.append(f"[{label}]")
        for key, value in summary.items():
            lines.append(f"{key}: {value}")
        lines.append("")
    summary_path.write_text("\\n".join(lines))

    for label, summary in summaries.items():
        print(label, summary)


if __name__ == "__main__":
    main()
