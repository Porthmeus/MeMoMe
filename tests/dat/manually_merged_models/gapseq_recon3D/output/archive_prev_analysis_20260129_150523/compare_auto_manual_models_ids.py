#!/usr/bin/env python3
"""Compare auto vs manual merged models after canonicalizing IDs."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import cobra
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[5]
sys.path.append(str(REPO_ROOT))
from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix  # noqa: E402

SPECIAL_PREFIXES = ("EX_", "DM_", "SK_", "sink_")


def canonical_met_id(met_id: str) -> str:
    base = met_id
    for prefix in ("model1_", "model2_", "H_", "M_"):
        if base.startswith(prefix):
            base = base[len(prefix) :]
            break
    cleaned = handle_metabolites_prefix_suffix(base)
    return cleaned if cleaned is not None else base


def strip_reaction_compartment(rxn_id: str) -> str:
    rid = rxn_id
    if rid.endswith("_t"):
        rid = rid[:-2]
    rid = re.sub(r"\([A-Za-z0-9_]+\)$", "", rid)
    rid = re.sub(r"\[[A-Za-z0-9_]+\]$", "", rid)
    rid = re.sub(r"_[a-z]\d?$", "", rid)
    return rid


def _strip_model_prefix(rid: str) -> str:
    for prefix in ("model1_", "model2_", "H_", "M_"):
        if rid.startswith(prefix):
            return rid[len(prefix) :]
    return rid


def _canonical_base_id(base: str) -> str:
    base = strip_reaction_compartment(base)
    cleaned = handle_metabolites_prefix_suffix(base)
    return cleaned if cleaned is not None else strip_reaction_compartment(base)


def canonical_reaction_id(rxn_id: str) -> str:
    rid = rxn_id

    if rid.startswith("TR_model1_") or rid.startswith("TR_model2_"):
        rid = "TR_" + rid.split("_", 2)[2]

    if rid.startswith("TR_"):
        base = _strip_model_prefix(rid[len("TR_") :])
        return "TR_" + _canonical_base_id(base)
    if rid.startswith(("IEX_", "H_IEX_", "M_IEX_")):
        if rid.startswith("IEX_"):
            base = _strip_model_prefix(rid[len("IEX_") :])
        else:
            base = _strip_model_prefix(rid.split("IEX_", 1)[1])
        return "TR_" + _canonical_base_id(base)

    stripped = _strip_model_prefix(rid)
    if stripped.startswith("IEX_"):
        return "TR_" + _canonical_base_id(stripped[len("IEX_") :])

    for sp in SPECIAL_PREFIXES:
        if rid.startswith(sp):
            base = _strip_model_prefix(rid[len(sp) :])
            return sp + _canonical_base_id(base)

    # Exchange reactions: normalize EX_model*_X_t and H_EX_*[f] to EX_<base>.
    ex_prefixes = ("EX_", "H_EX_", "M_EX_")
    for ex_prefix in ex_prefixes:
        if rid.startswith(ex_prefix):
            base = rid[len(ex_prefix) :]
            base = _strip_model_prefix(base)
            return "EX_" + _canonical_base_id(base)
    if rid.startswith("EX_model1_") or rid.startswith("EX_model2_"):
        base = rid.split("_", 2)[2]
        return "EX_" + _canonical_base_id(base)

    rid = _strip_model_prefix(rid)
    for sp in SPECIAL_PREFIXES:
        if rid.startswith(sp):
            return sp + _canonical_base_id(rid[len(sp) :])
    return _canonical_base_id(rid)


def exchange_reaction_ids(model: cobra.Model) -> list[cobra.Reaction]:
    exchange = []
    for rxn in model.reactions:
        if rxn.id.startswith("EX_") or rxn.id.startswith("H_EX_") or rxn.id.startswith("M_EX_"):
            exchange.append(rxn)
    return exchange


def exchange_bounds_map(model: cobra.Model) -> dict[str, list[dict[str, str | float]]]:
    bounds: dict[str, list[dict[str, str | float]]] = {}
    for rxn in exchange_reaction_ids(model):
        if not rxn.metabolites:
            continue
        met = next(iter(rxn.metabolites))
        base = canonical_met_id(met.id)
        bounds.setdefault(base, []).append(
            {
                "rxn_id": rxn.id,
                "rxn_name": rxn.name or "",
                "lb": rxn.lower_bound,
                "met_id": met.id,
            }
        )
    return bounds


def ex_t_bounds_map(model: cobra.Model) -> dict[str, list[dict[str, str | float]]]:
    bounds: dict[str, list[dict[str, str | float]]] = {}
    for rxn in model.reactions:
        if not (rxn.id.startswith("EX_") and rxn.id.endswith("_t")):
            continue
        if not rxn.metabolites:
            continue
        met = next(iter(rxn.metabolites))
        base = canonical_met_id(met.id)
        bounds.setdefault(base, []).append(
            {
                "rxn_id": rxn.id,
                "rxn_name": rxn.name or "",
                "lb": rxn.lower_bound,
                "met_id": met.id,
            }
        )
    return bounds


def write_set(output_dir: Path, name: str, values: set[str]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(sorted(values), columns=[name])
    df.to_csv(output_dir / f"{name}.csv", index=False)


def write_diff_sets(output_dir: Path, prefix: str, auto_set: set[str], manual_set: set[str]) -> None:
    write_set(output_dir, f"{prefix}_intersection", auto_set & manual_set)
    write_set(output_dir, f"{prefix}_only_auto", auto_set - manual_set)
    write_set(output_dir, f"{prefix}_only_manual", manual_set - auto_set)


def compare_diets(
    output_dir: Path,
    auto_bounds: dict[str, list[dict[str, str | float]]],
    manual_bounds: dict[str, list[dict[str, str | float]]],
    tag: str,
) -> None:
    auto_keys = set(auto_bounds)
    manual_keys = set(manual_bounds)

    overlap_rows = []
    for key in sorted(auto_keys & manual_keys):
        auto_vals = auto_bounds[key]
        manual_vals = manual_bounds[key]
        overlap_rows.append(
            {
                "metabolite": key,
                "auto_ids": ";".join(str(r["rxn_id"]) for r in auto_vals),
                "auto_names": ";".join(str(r["rxn_name"]) for r in auto_vals),
                "auto_lbs": ";".join(str(r["lb"]) for r in auto_vals),
                "manual_ids": ";".join(str(r["rxn_id"]) for r in manual_vals),
                "manual_names": ";".join(str(r["rxn_name"]) for r in manual_vals),
                "manual_lbs": ";".join(str(r["lb"]) for r in manual_vals),
            }
        )
    pd.DataFrame(overlap_rows).to_csv(output_dir / f"diet_{tag}_overlap.csv", index=False)

    write_set(output_dir, f"diet_{tag}_only_auto", auto_keys - manual_keys)
    write_set(output_dir, f"diet_{tag}_only_manual", manual_keys - auto_keys)


def tr_shuttled_map(model: cobra.Model) -> dict[str, list[dict[str, str]]]:
    shuttled: dict[str, list[dict[str, str]]] = {}
    for rxn in model.reactions:
        if not rxn.id.startswith("TR_"):
            continue
        translation = [met for met in rxn.metabolites if met.compartment == "t" or met.id.endswith("_t")]
        if translation:
            base = canonical_met_id(translation[0].id)
            shuttled.setdefault(base, []).append(
                {"rxn_id": rxn.id, "rxn_name": rxn.name or ""}
            )
    return shuttled


def iex_shuttled_map(model: cobra.Model) -> dict[str, list[dict[str, str]]]:
    shuttled: dict[str, list[dict[str, str]]] = {}
    for rxn in model.reactions:
        if "IEX_" not in rxn.id:
            continue
        feed = [met for met in rxn.metabolites if met.compartment == "f" or met.id.endswith("[f]")]
        if feed:
            base = canonical_met_id(feed[0].id)
            shuttled.setdefault(base, []).append(
                {"rxn_id": rxn.id, "rxn_name": rxn.name or ""}
            )
    return shuttled


def _build_met_info(model: cobra.Model) -> dict[str, list[dict[str, str]]]:
    info: dict[str, list[dict[str, str]]] = {}
    for met in model.metabolites:
        base = canonical_met_id(met.id)
        info.setdefault(base, []).append(
            {
                "met_id": met.id,
                "name": met.name or "",
                "compartment": met.compartment or "",
            }
        )
    return info


def _build_rxn_info(model: cobra.Model) -> dict[str, list[dict[str, str]]]:
    info: dict[str, list[dict[str, str]]] = {}
    for rxn in model.reactions:
        base = canonical_reaction_id(rxn.id)
        info.setdefault(base, []).append({"rxn_id": rxn.id, "name": rxn.name or ""})
    return info


def _write_only_with_names(
    output_dir: Path,
    filename: str,
    only_set: set[str],
    info: dict[str, list[dict[str, str]]],
    kind: str,
) -> None:
    rows = []
    for key in sorted(only_set):
        entries = info.get(key, [])
        rows.append(
            {
                "canonical_id": key,
                "ids": ";".join(e["met_id" if kind == "met" else "rxn_id"] for e in entries),
                "names": ";".join(e["name"] for e in entries),
                "compartments": ";".join(e.get("compartment", "") for e in entries) if kind == "met" else "",
            }
        )
    pd.DataFrame(rows).to_csv(output_dir / filename, index=False)


def _write_diet_diff(
    output_dir: Path,
    auto_bounds: dict[str, list[dict[str, str | float]]],
    manual_bounds: dict[str, list[dict[str, str | float]]],
) -> None:
    rows = []
    for key in sorted(set(auto_bounds) & set(manual_bounds)):
        auto_vals = auto_bounds[key]
        manual_vals = manual_bounds[key]
        auto_lbs = [float(v["lb"]) for v in auto_vals]
        manual_lbs = [float(v["lb"]) for v in manual_vals]
        if sorted(auto_lbs) == sorted(manual_lbs):
            continue
        rows.append(
            {
                "metabolite": key,
                "auto_ids": ";".join(str(r["rxn_id"]) for r in auto_vals),
                "auto_names": ";".join(str(r["rxn_name"]) for r in auto_vals),
                "auto_lbs": ";".join(str(r["lb"]) for r in auto_vals),
                "manual_ids": ";".join(str(r["rxn_id"]) for r in manual_vals),
                "manual_names": ";".join(str(r["rxn_name"]) for r in manual_vals),
                "manual_lbs": ";".join(str(r["lb"]) for r in manual_vals),
            }
        )
    pd.DataFrame(rows).to_csv(output_dir / "diet_lb_differences.csv", index=False)


def _write_shuttled_only(
    output_dir: Path,
    filename: str,
    only_set: set[str],
    shuttled_map: dict[str, list[dict[str, str]]],
) -> None:
    rows = []
    for key in sorted(only_set):
        entries = shuttled_map.get(key, [])
        rows.append(
            {
                "metabolite": key,
                "rxn_ids": ";".join(e["rxn_id"] for e in entries),
                "rxn_names": ";".join(e["rxn_name"] for e in entries),
            }
        )
    pd.DataFrame(rows).to_csv(output_dir / filename, index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare auto vs manual merged model IDs.")
    parser.add_argument(
        "--auto-path",
        type=Path,
        default=REPO_ROOT
        / "tests/dat/manually_merged_models/gapseq_recon3D/output/automatically_merged_metamodel.xml",
    )
    parser.add_argument(
        "--manual-path",
        type=Path,
        default=REPO_ROOT
        / "tests/dat/manually_merged_models/gapseq_recon3D/output/manual_diet/merged_model_2025_prefixed_normalized_diet.xml",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "compare_models_ids",
    )
    args = parser.parse_args()

    auto_model = cobra.io.read_sbml_model(str(args.auto_path))
    manual_model = cobra.io.read_sbml_model(str(args.manual_path))

    auto_mets = {canonical_met_id(m.id) for m in auto_model.metabolites}
    manual_mets = {canonical_met_id(m.id) for m in manual_model.metabolites}
    write_diff_sets(args.output_dir, "metabolites", auto_mets, manual_mets)

    auto_rxns = {canonical_reaction_id(r.id) for r in auto_model.reactions}
    manual_rxns = {canonical_reaction_id(r.id) for r in manual_model.reactions}
    write_diff_sets(args.output_dir, "reactions", auto_rxns, manual_rxns)

    # Diet bounds: auto uses EX_*_t, manual uses H_EX_* (and other EX variants).
    auto_diet = ex_t_bounds_map(auto_model)
    manual_diet = exchange_bounds_map(manual_model)
    compare_diets(args.output_dir, auto_diet, manual_diet, "auto_ex_t_vs_manual_ex")

    # TR vs IEX shuttled metabolites.
    auto_tr_map = tr_shuttled_map(auto_model)
    manual_iex_map = iex_shuttled_map(manual_model)
    auto_tr_mets = set(auto_tr_map)
    manual_iex_mets = set(manual_iex_map)
    write_diff_sets(args.output_dir, "shuttled", auto_tr_mets, manual_iex_mets)

    summary_dir = args.output_dir / "summary_diffs"
    summary_dir.mkdir(parents=True, exist_ok=True)

    auto_met_info = _build_met_info(auto_model)
    manual_met_info = _build_met_info(manual_model)
    _write_only_with_names(
        summary_dir,
        "metabolites_only_auto_with_names.csv",
        auto_mets - manual_mets,
        auto_met_info,
        "met",
    )
    _write_only_with_names(
        summary_dir,
        "metabolites_only_manual_with_names.csv",
        manual_mets - auto_mets,
        manual_met_info,
        "met",
    )

    auto_rxn_info = _build_rxn_info(auto_model)
    manual_rxn_info = _build_rxn_info(manual_model)
    _write_only_with_names(
        summary_dir,
        "reactions_only_auto_with_names.csv",
        auto_rxns - manual_rxns,
        auto_rxn_info,
        "rxn",
    )
    _write_only_with_names(
        summary_dir,
        "reactions_only_manual_with_names.csv",
        manual_rxns - auto_rxns,
        manual_rxn_info,
        "rxn",
    )

    _write_only_with_names(
        summary_dir,
        "diet_only_auto_with_names.csv",
        set(auto_diet) - set(manual_diet),
        auto_met_info,
        "met",
    )
    _write_only_with_names(
        summary_dir,
        "diet_only_manual_with_names.csv",
        set(manual_diet) - set(auto_diet),
        manual_met_info,
        "met",
    )

    _write_diet_diff(summary_dir, auto_diet, manual_diet)

    _write_shuttled_only(
        summary_dir,
        "shuttled_only_auto_with_names.csv",
        auto_tr_mets - manual_iex_mets,
        auto_tr_map,
    )
    _write_shuttled_only(
        summary_dir,
        "shuttled_only_manual_with_names.csv",
        manual_iex_mets - auto_tr_mets,
        manual_iex_map,
    )

    summary = [
        f"auto metabolites: {len(auto_mets)}",
        f"manual metabolites: {len(manual_mets)}",
        f"metabolite overlap: {len(auto_mets & manual_mets)}",
        f"auto reactions: {len(auto_rxns)}",
        f"manual reactions: {len(manual_rxns)}",
        f"reaction overlap: {len(auto_rxns & manual_rxns)}",
        f"auto diet entries: {len(auto_diet)}",
        f"manual diet entries: {len(manual_diet)}",
        f"diet overlap: {len(set(auto_diet) & set(manual_diet))}",
        f"auto TR shuttled: {len(auto_tr_mets)}",
        f"manual IEX shuttled: {len(manual_iex_mets)}",
        f"shuttled overlap: {len(auto_tr_mets & manual_iex_mets)}",
    ]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "summary.txt").write_text("\n".join(summary) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
