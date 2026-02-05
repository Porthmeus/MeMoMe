#!/usr/bin/env python3
"""Apply auto-model EX_t bounds to manual EX reactions and compare diets."""

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
    for suffix in ("_e", "_c", "_p", "_b", "_f", "_e0", "_c0", "_p0"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    return base


def _reaction_metabolite(rxn: cobra.Reaction) -> cobra.Metabolite | None:
    if not rxn.metabolites:
        return None
    return next(iter(rxn.metabolites))


def _auto_ex_t_bounds(auto_model_path: Path) -> dict[str, tuple[float, float]]:
    auto_model = cobra.io.read_sbml_model(str(auto_model_path))
    bounds: dict[str, tuple[float, float]] = {}
    for rxn in auto_model.reactions:
        if not (rxn.id.startswith("EX_") and rxn.id.endswith("_t")):
            continue
        met = _reaction_metabolite(rxn)
        if met is None:
            continue
        base = _canonical_met_id(met.id)
        bounds[base] = (rxn.lower_bound, rxn.upper_bound)
    return bounds


def _manual_ex_candidates(model: cobra.Model) -> list[cobra.Reaction]:
    candidates = []
    for rxn in model.reactions:
        if not rxn.id.startswith(("H_EX_", "M_EX_", "EX_")):
            continue
        candidates.append(rxn)
    return candidates


def apply_manual_diet(
    model: cobra.Model, auto_model_path: Path
) -> tuple[int, int]:
    auto_bounds = _auto_ex_t_bounds(auto_model_path)
    candidates = 0
    applied = 0
    for rxn in _manual_ex_candidates(model):
        candidates += 1
        met = _reaction_metabolite(rxn)
        if met is None:
            continue
        base = _canonical_met_id(met.id)
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
    parser.add_argument(
        "--diet-compare-dir",
        type=Path,
        default=Path(
            "tests/dat/manually_merged_models/gapseq_recon3D/output/manual_diet/"
            "diet_compare"
        ),
        help="Output directory for diet comparison CSVs.",
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

    # Compare diets between auto EX_t and manual EX reactions.
    auto_bounds = _auto_ex_t_bounds(args.auto_model)
    manual_bounds: dict[str, tuple[float, float]] = {}
    for rxn in _manual_ex_candidates(model):
        met = _reaction_metabolite(rxn)
        if met is None:
            continue
        base = _canonical_met_id(met.id)
        manual_bounds[base] = (rxn.lower_bound, rxn.upper_bound)

    auto_set = set(auto_bounds)
    manual_set = set(manual_bounds)
    overlap = sorted(auto_set & manual_set)
    only_auto = sorted(auto_set - manual_set)
    only_manual = sorted(manual_set - auto_set)

    args.diet_compare_dir.mkdir(parents=True, exist_ok=True)

    overlap_lines = ["met_id,auto_lb,auto_ub,manual_lb,manual_ub,lb_diff,ub_diff"]
    for met in overlap:
        auto_lb, auto_ub = auto_bounds[met]
        manual_lb, manual_ub = manual_bounds[met]
        overlap_lines.append(
            f"{met},{auto_lb},{auto_ub},{manual_lb},{manual_ub},"
            f"{manual_lb - auto_lb},{manual_ub - auto_ub}"
        )
    (args.diet_compare_dir / "diet_overlap.csv").write_text("\n".join(overlap_lines) + "\n")
    only_auto_lines = ["met_id,auto_lb,auto_ub"]
    for met in only_auto:
        auto_lb, auto_ub = auto_bounds[met]
        only_auto_lines.append(f"{met},{auto_lb},{auto_ub}")
    (args.diet_compare_dir / "diet_only_auto.csv").write_text("\n".join(only_auto_lines) + "\n")

    only_manual_lines = ["met_id,manual_lb,manual_ub"]
    for met in only_manual:
        manual_lb, manual_ub = manual_bounds[met]
        only_manual_lines.append(f"{met},{manual_lb},{manual_ub}")
    (args.diet_compare_dir / "diet_only_manual.csv").write_text("\n".join(only_manual_lines) + "\n")


if __name__ == "__main__":
    main()
