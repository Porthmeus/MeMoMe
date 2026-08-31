#!/usr/bin/env python3
"""Apply auto-model EX_t bounds to manual H_EX reactions and save a new SBML."""

from __future__ import annotations

import argparse
from pathlib import Path

import cobra


def _canonical_met_id(met_id: str) -> str:
    base = met_id
    for prefix in ("H_", "M_", "model1_", "model2_"):
        if base.startswith(prefix):
            base = base[len(prefix) :]
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
    bounds: dict[str, tuple[float, float]] = {}
    for rxn in auto_model.reactions:
        if not (rxn.id.startswith("EX_") and rxn.id.endswith("_t")):
            continue
        if not rxn.metabolites:
            continue
        met = next(iter(rxn.metabolites))
        base = _canonical_met_id(met.id)
        bounds[base] = (rxn.lower_bound, rxn.upper_bound)
    return bounds


def apply_manual_diet(
    model: cobra.Model, auto_model_path: Path
) -> tuple[int, int]:
    auto_bounds = _auto_ex_t_bounds(auto_model_path)
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
    return applied, candidates


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Apply auto-model EX_t bounds to manual H_EX reactions."
    )
    parser.add_argument(
        "--auto-model",
        type=Path,
        default=Path(
            "tests/dat/manually_merged_models/gapseq_recon3D/output/automatically_merged_metamodel.xml"
        ),
        help="Auto model used to transfer EX_t diet to manual H_EX reactions.",
    )
    parser.add_argument(
        "--manual-input",
        type=Path,
        default=Path(
            "tests/dat/manually_merged_models/gapseq_recon3D/merged_model_2025_prefixed_normalized.xml"
        ),
        help="Manual model to update with auto-derived diet bounds.",
    )
    parser.add_argument(
        "--manual-output",
        type=Path,
        default=Path(
            "tests/dat/manually_merged_models/gapseq_recon3D/output/manual_diet/"
            "merged_model_2025_prefixed_normalized_diet.xml"
        ),
        help="Output SBML with diet bounds applied.",
    )
    parser.add_argument(
        "--report-path",
        type=Path,
        default=Path(
            "tests/dat/manually_merged_models/gapseq_recon3D/output/manual_diet/"
            "manual_diet_applied.txt"
        ),
        help="Text report with diet application counts.",
    )
    args = parser.parse_args()

    model = cobra.io.read_sbml_model(str(args.manual_input))
    applied, candidates = apply_manual_diet(model, args.auto_model)

    args.manual_output.parent.mkdir(parents=True, exist_ok=True)
    cobra.io.write_sbml_model(model, str(args.manual_output))

    args.report_path.parent.mkdir(parents=True, exist_ok=True)
    args.report_path.write_text(
        f"applied {applied} of {candidates} H_EX bounds from {args.auto_model}\n"
    )

    print(
        f"Saved dieted manual model to {args.manual_output} "
        f"(applied {applied} of {candidates} candidates)."
    )


if __name__ == "__main__":
    main()
