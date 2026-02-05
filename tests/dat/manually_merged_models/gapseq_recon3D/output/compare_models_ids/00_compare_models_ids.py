#!/usr/bin/env python3
"""Compare auto vs manual merged models after canonicalizing IDs."""

from __future__ import annotations

import re
from pathlib import Path

import cobra
import pandas as pd

from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix


SPECIAL_PREFIXES = ("EX_", "DM_", "SK_", "sink_")

COMPARTMENT_SUFFIXES = [
    "(e)", "(c)", "(p)", "(m)", "(x)", "(n)", "(r)", "(g)", "(l)", "(u)", "(b)", "(f)", "(t)",
    "[e]", "[c]", "[p]", "[m]", "[x]", "[n]", "[r]", "[g]", "[l]", "[u]", "[b]", "[f]", "[t]",
    "_e", "_c", "_p", "_m", "_x", "_n", "_r", "_g", "_l", "_u", "_b", "_f", "_t",
    "_e0", "_c0", "_p0",
]


def strip_compartment_suffix(value: str) -> str:
    for suffix in COMPARTMENT_SUFFIXES:
        if value.endswith(suffix):
            return value[: -len(suffix)]
    return value


def strip_reaction_compartment(rxn_id: str) -> str:
    rid = rxn_id
    if rid.endswith("_t"):
        rid = rid[:-2]
    rid = re.sub(r"\([A-Za-z0-9_]+\)$", "", rid)
    rid = re.sub(r"\[[A-Za-z0-9_]+\]$", "", rid)
    rid = re.sub(r"_[a-z]\d?$", "", rid)
    return rid


def strip_model_prefix(value: str) -> str:
    for prefix in ("model1_", "model2_", "H_", "M_"):
        if value.startswith(prefix):
            return value[len(prefix) :]
    return value


def canonical_met_id(met_id: str) -> str:
    base = strip_model_prefix(met_id)
    cleaned = handle_metabolites_prefix_suffix(base)
    return cleaned if cleaned is not None else base


def canonical_reaction_id(rxn_id: str) -> str:
    rid = rxn_id

    if rid.startswith("TR_model1_") or rid.startswith("TR_model2_"):
        rid = "TR_" + rid.split("_", 2)[2]

    if rid.startswith("TR_"):
        base = strip_model_prefix(rid[len("TR_") :])
        return "TR_" + strip_reaction_compartment(base)

    if rid.startswith(("IEX_", "H_IEX_", "M_IEX_")):
        if rid.startswith("IEX_"):
            base = strip_model_prefix(rid[len("IEX_") :])
        else:
            base = strip_model_prefix(rid.split("IEX_", 1)[1])
        return "TR_" + strip_reaction_compartment(base)

    stripped = strip_model_prefix(rid)
    if stripped.startswith("IEX_"):
        return "TR_" + strip_reaction_compartment(stripped[len("IEX_") :])

    for prefix in SPECIAL_PREFIXES:
        if rid.startswith(prefix):
            base = strip_model_prefix(rid[len(prefix) :])
            return prefix + strip_reaction_compartment(base)

    ex_prefixes = ("EX_", "H_EX_", "M_EX_")
    for ex_prefix in ex_prefixes:
        if rid.startswith(ex_prefix):
            base = rid[len(ex_prefix) :]
            base = strip_model_prefix(base)
            return "EX_" + strip_reaction_compartment(base)

    if rid.startswith("EX_model1_") or rid.startswith("EX_model2_"):
        base = rid.split("_", 2)[2]
        return "EX_" + strip_reaction_compartment(base)

    rid = strip_model_prefix(rid)
    for prefix in SPECIAL_PREFIXES:
        if rid.startswith(prefix):
            return prefix + strip_reaction_compartment(rid[len(prefix) :])

    return strip_reaction_compartment(rid)


def exchange_reaction_ids(model: cobra.Model) -> list[cobra.Reaction]:
    exchange = []
    for rxn in model.reactions:
        if rxn.id.startswith("EX_") or rxn.id.startswith("H_EX_") or rxn.id.startswith("M_EX_"):
            exchange.append(rxn)
    return exchange


def write_set(output_dir: Path, name: str, values: set[str]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(sorted(values), columns=[name])
    df.to_csv(output_dir / f"{name}.csv", index=False)


def write_diff_sets(output_dir: Path, prefix: str, auto_set: set[str], manual_set: set[str]) -> None:
    write_set(output_dir, f"{prefix}_intersection", auto_set & manual_set)
    write_set(output_dir, f"{prefix}_only_auto", auto_set - manual_set)
    write_set(output_dir, f"{prefix}_only_manual", manual_set - auto_set)


def main() -> None:
    output_dir = Path(__file__).resolve().parent
    auto_model_path = output_dir.parent / "automatically_merged_metamodel.xml"
    manual_model_path = output_dir.parent / "merged_model_2025_prefixed_normalized_diet_fixIEX.xml"

    auto_model = cobra.io.read_sbml_model(str(auto_model_path))
    manual_model = cobra.io.read_sbml_model(str(manual_model_path))

    auto_mets = {canonical_met_id(m.id) for m in auto_model.metabolites}
    manual_mets = {canonical_met_id(m.id) for m in manual_model.metabolites}
    write_diff_sets(output_dir, "metabolites", auto_mets, manual_mets)

    auto_rxns = {canonical_reaction_id(r.id) for r in auto_model.reactions}
    manual_rxns = {canonical_reaction_id(r.id) for r in manual_model.reactions}
    write_diff_sets(output_dir, "reactions", auto_rxns, manual_rxns)

    auto_diet_bounds = {}
    for r in auto_model.reactions:
        if r.id.startswith("EX_") and r.id.endswith("_t") and r.metabolites:
            met_id = canonical_met_id(next(iter(r.metabolites)).id)
            auto_diet_bounds[met_id] = (r.lower_bound, r.upper_bound)

    manual_diet_bounds = {}
    for r in exchange_reaction_ids(manual_model):
        if r.metabolites:
            met_id = canonical_met_id(next(iter(r.metabolites)).id)
            manual_diet_bounds[met_id] = (r.lower_bound, r.upper_bound)

    auto_diet = set(auto_diet_bounds)
    manual_diet = set(manual_diet_bounds)
    write_diff_sets(output_dir, "diet", auto_diet, manual_diet)

    # Diet intersection with bounds.
    overlap = sorted(auto_diet & manual_diet)
    overlap_lines = ["met_id,auto_lb,auto_ub,manual_lb,manual_ub,lb_diff,ub_diff"]
    for met in overlap:
        auto_lb, auto_ub = auto_diet_bounds[met]
        manual_lb, manual_ub = manual_diet_bounds[met]
        overlap_lines.append(
            f"{met},{auto_lb},{auto_ub},{manual_lb},{manual_ub},"
            f"{manual_lb - auto_lb},{manual_ub - auto_ub}"
        )
    (output_dir / "diet_intersection.csv").write_text("\n".join(overlap_lines) + "\n")

    # Diet only auto with bounds.
    only_auto = sorted(auto_diet - manual_diet)
    only_auto_lines = ["met_id,auto_lb,auto_ub"]
    for met in only_auto:
        auto_lb, auto_ub = auto_diet_bounds[met]
        only_auto_lines.append(f"{met},{auto_lb},{auto_ub}")
    (output_dir / "diet_only_auto.csv").write_text("\n".join(only_auto_lines) + "\n")

    # Diet only manual with bounds.
    only_manual = sorted(manual_diet - auto_diet)
    only_manual_lines = ["met_id,manual_lb,manual_ub"]
    for met in only_manual:
        manual_lb, manual_ub = manual_diet_bounds[met]
        only_manual_lines.append(f"{met},{manual_lb},{manual_ub}")
    (output_dir / "diet_only_manual.csv").write_text("\n".join(only_manual_lines) + "\n")

    auto_tr_mets = {
        canonical_met_id(m.id)
        for r in auto_model.reactions
        if r.id.startswith("TR_")
        for m in r.metabolites
        if m.compartment == "t" or m.id.endswith("_t")
    }
    manual_iex_mets = {
        canonical_met_id(m.id)
        for r in manual_model.reactions
        if "IEX_" in r.id
        for m in r.metabolites
        if m.compartment == "f" or m.id.endswith("[f]")
    }
    write_diff_sets(output_dir, "shuttled", auto_tr_mets, manual_iex_mets)


if __name__ == "__main__":
    main()
