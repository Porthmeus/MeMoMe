#!/usr/bin/env python3
"""Run FVA computations for auto/manual models (no postprocessing)."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import os
from pathlib import Path
from typing import Iterable

import cobra
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


def _is_host_reaction_manual(rxn: cobra.Reaction) -> bool:
    if not rxn.id.startswith("H_"):
        return False
    if rxn.id.startswith(("H_EX_", "H_IEX_")):
        return False
    if any(met.compartment == "f" for met in rxn.metabolites):
        return False
    return True


def _host_reaction_ids(model: cobra.Model, spec: ModelSpec, max_reactions: int | None) -> list[str]:
    if spec.label == "auto":
        host_ids = [rxn.id for rxn in model.reactions if _is_host_reaction_auto(rxn.id)]
    else:
        host_ids = [rxn.id for rxn in model.reactions if _is_host_reaction_manual(rxn)]
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
    # Force zero exchange of O2 for microbe in both auto and manual models.
    if spec.label == "auto":
        return [rxn.id for rxn in model.reactions if "TR_model2_cpd00007" in rxn.id]
    forced = []
    for rxn in model.reactions:
        if "IEX_" not in rxn.id or "cpd00007" not in rxn.id:
            continue
        if any(met.id.startswith("M_") for met in rxn.metabolites):
            forced.append(rxn.id)
    return forced


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


def run_model_fvas(
    spec: ModelSpec,
    output_dir: Path,
    processes: int,
    max_reactions: int | None,
    max_exchanges: int | None,
) -> dict[str, float | int]:
    model = cobra.io.read_sbml_model(str(spec.path))
    host_reaction_ids = _host_reaction_ids(model, spec, max_reactions)
    biomass_ids = _identify_biomass_ids(model, spec.host_prefix, spec.micro_prefix)
    original_lbs = {
        key: model.reactions.get_by_id(rxn_id).lower_bound
        for key, rxn_id in biomass_ids.items()
    }
    forced_zero_ids = _forced_zero_exchange_ids(model, spec)

    output_dir.mkdir(parents=True, exist_ok=True)

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
    try:
        baseline_fva = _run_fva(baseline_model, host_reaction_ids, processes)
    except Exception as exc:  # noqa: BLE001
        raise Infeasible(f"{spec.label} baseline model infeasible") from exc
    baseline_fva.to_csv(output_dir / f"{spec.label}_host_fva_baseline.csv")

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
    if max_exchanges is not None:
        bacterial_exchange_ids = bacterial_exchange_ids[:max_exchanges]

    biomass_id = biomass_ids["host"]
    exchange_fva_ids = list(bacterial_exchange_ids)
    if biomass_id not in exchange_fva_ids:
        exchange_fva_ids.append(biomass_id)
    exchange_fva = _run_fva(host_only_model, exchange_fva_ids, processes)
    exchange_fva.to_csv(output_dir / f"{spec.label}_bacterial_exchange_fva.csv")

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
    _close_reactions(closed_model, bacterial_exchange_ids)
    try:
        closed_fva = _run_fva(closed_model, host_reaction_ids, processes)
    except Exception:  # noqa: BLE001
        pd.DataFrame({"status": ["infeasible"]}).to_csv(
            output_dir / f"{spec.label}_host_fva_closed.csv", index=False
        )
        run_info = {
            "host_reactions": len(host_reaction_ids),
            "bacterial_exchange_closed": len(bacterial_exchange_ids),
            "closed_feasible": False,
            "forced_zero_exchange_count": len(forced_zero_ids),
        }
        print(spec.label, run_info)
        return run_info

    closed_fva.to_csv(output_dir / f"{spec.label}_host_fva_closed.csv")
    run_info = {
        "host_reactions": len(host_reaction_ids),
        "bacterial_exchange_closed": len(bacterial_exchange_ids),
        "closed_feasible": True,
        "forced_zero_exchange_count": len(forced_zero_ids),
    }
    print(spec.label, run_info)
    return run_info


def main() -> None:
    parser = argparse.ArgumentParser(description="Run FVA computations for auto/manual models.")
    parser.add_argument(
        "--auto-path",
        type=Path,
        default=Path(
            "tests/dat/manually_merged_models/gapseq_recon3D/output/automatically_merged_metamodel.xml"
        ),
    )
    parser.add_argument(
        "--manual-path",
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
    parser.add_argument(
        "--max-reactions",
        type=int,
        default=None,
        help="Optional cap on host reactions (debugging only).",
    )
    parser.add_argument(
        "--max-exchanges",
        type=int,
        default=None,
        help="Optional cap on bacterial exchange reactions (debugging only).",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=None,
        help="Number of worker processes for FVA (defaults to CPU count).",
    )
    args = parser.parse_args()

    if args.processes is None:
        processes = max(1, os.cpu_count() or 1)
    else:
        processes = max(1, args.processes)

    specs = [
        ModelSpec("auto", args.auto_path, "model1_", "model2_"),
        ModelSpec("manual", args.manual_path, "H_", "M_"),
    ]
    for spec in specs:
        run_model_fvas(
            spec,
            args.output_dir,
            processes,
            args.max_reactions,
            args.max_exchanges,
        )


if __name__ == "__main__":
    main()
