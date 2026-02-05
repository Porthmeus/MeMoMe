#!/usr/bin/env python3
"""Detect reaction ID mismatches caused by naming conventions.

Produces:
- reaction_only_same_metabolite.csv
- residual_mismatches_with_names.csv
- ex_iex_feed_mismatch_refined.csv
"""
from __future__ import annotations

from pathlib import Path
import importlib.util
import csv
import cobra

BASE = Path(__file__).resolve().parent
AUTO_MODEL = BASE.parent / "automatically_merged_metamodel.xml"
MANUAL_MODEL = BASE.parent / "merged_model_2025_prefixed_normalized_diet_fixIEX.xml"


def _load_canonicalizers():
    src = BASE / "00_compare_models_ids.py"
    spec = importlib.util.spec_from_file_location("cmp", src)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.canonical_reaction_id, mod.canonical_met_id


def _build_rxn_map(model, canon_rxn):
    mapping = {}
    for rxn in model.reactions:
        key = canon_rxn(rxn.id)
        if key not in mapping:
            mapping[key] = rxn
    return mapping


def _rxn_met_canon(rxn, canon_met):
    if rxn is None:
        return None
    met = next(iter(rxn.metabolites), None)
    return canon_met(met.id) if met is not None else None


def _load_list(path: Path):
    return [row.strip() for row in path.read_text().splitlines()[1:] if row.strip()]


def _met_in_comp(met, comp: str) -> bool:
    if getattr(met, "compartment", None) == comp:
        return True
    return met.id.endswith(f"[{comp}]") or met.id.endswith(f"_{comp}")


def main() -> None:
    canon_rxn, canon_met = _load_canonicalizers()

    auto = cobra.io.read_sbml_model(str(AUTO_MODEL))
    manual = cobra.io.read_sbml_model(str(MANUAL_MODEL))

    auto_only = _load_list(BASE / "reactions_only_auto.csv")
    manual_only = _load_list(BASE / "reactions_only_manual.csv")

    auto_map = _build_rxn_map(auto, canon_rxn)
    manual_map = _build_rxn_map(manual, canon_rxn)

    auto_only_by_met = {}
    for rid in auto_only:
        cm = _rxn_met_canon(auto_map.get(rid), canon_met)
        if cm is None:
            continue
        auto_only_by_met.setdefault(cm, []).append(rid)

    manual_only_by_met = {}
    for rid in manual_only:
        cm = _rxn_met_canon(manual_map.get(rid), canon_met)
        if cm is None:
            continue
        manual_only_by_met.setdefault(cm, []).append(rid)

    shared_mets = sorted(set(auto_only_by_met) & set(manual_only_by_met))

    same_met = BASE / "reaction_only_same_metabolite.csv"
    lines = ["met_id,auto_reaction_ids,manual_reaction_ids"]
    for met in shared_mets:
        lines.append(f"{met}," + "|".join(auto_only_by_met[met]) + "," + "|".join(manual_only_by_met[met]))
    same_met.write_text("\n".join(lines) + "\n")

    name_mismatch = set()
    for met in shared_mets:
        name_mismatch.update(auto_only_by_met[met])
        name_mismatch.update(manual_only_by_met[met])

    residual_auto = [rid for rid in auto_only if rid not in name_mismatch]
    residual_manual = [rid for rid in manual_only if rid not in name_mismatch]

    residual_out = BASE / "residual_mismatches_with_names.csv"
    lines = ["model,canonical_reaction_id,actual_reaction_id,met_id,met_name"]
    for rid in residual_auto:
        rxn = auto_map.get(rid)
        met = next(iter(rxn.metabolites), None) if rxn else None
        lines.append(f"auto,{rid},{rxn.id if rxn else ''},{met.id if met else ''},{met.name if met else ''}")
    for rid in residual_manual:
        rxn = manual_map.get(rid)
        met = next(iter(rxn.metabolites), None) if rxn else None
        lines.append(f"manual,{rid},{rxn.id if rxn else ''},{met.id if met else ''},{met.name if met else ''}")
    residual_out.write_text("\n".join(lines) + "\n")

    # EX -> IEX/TR refined feed check (compartments)
    auto_ex = {}
    for r in auto.reactions:
        if r.id.startswith("EX_") and r.id.endswith("_t") and r.metabolites:
            m = next(iter(r.metabolites))
            if _met_in_comp(m, "t") or m.id.endswith("_t"):
                auto_ex.setdefault(canon_met(m.id), []).append(r.id)

    auto_tr = {}
    for r in auto.reactions:
        if r.id.startswith("TR_") and r.metabolites:
            for m in r.metabolites:
                if _met_in_comp(m, "t") or m.id.endswith("_t"):
                    auto_tr.setdefault(canon_met(m.id), []).append(r.id)

    manual_ex = {}
    for r in manual.reactions:
        if r.id.startswith("EX_") and r.metabolites:
            m = next(iter(r.metabolites))
            if _met_in_comp(m, "f"):
                manual_ex.setdefault(canon_met(m.id), []).append(r.id)

    manual_iex = {}
    for r in manual.reactions:
        if "IEX_" in r.id and r.metabolites:
            for m in r.metabolites:
                if _met_in_comp(m, "f"):
                    manual_iex.setdefault(canon_met(m.id), []).append(r.id)

    auto_ex_only = sorted(set(auto_ex) - set(auto_tr))
    manual_ex_only = sorted(set(manual_ex) - set(manual_iex))

    feed_out = BASE / "ex_iex_feed_mismatch_refined.csv"
    lines = ["model,met_id,ex_ids,link_ids"]
    for met in auto_ex_only:
        lines.append("auto," + met + "," + "|".join(auto_ex.get(met, [])) + "," + "|".join(auto_tr.get(met, [])))
    for met in manual_ex_only:
        lines.append("manual," + met + "," + "|".join(manual_ex.get(met, [])) + "," + "|".join(manual_iex.get(met, [])))
    feed_out.write_text("\n".join(lines) + "\n")

    print(f"same-met mismatches: {len(shared_mets)}")
    print(f"residual auto-only: {len(residual_auto)} | residual manual-only: {len(residual_manual)}")
    print(f"refined EX-TR/IEX mismatches: {len(auto_ex_only) + len(manual_ex_only)}")


if __name__ == "__main__":
    main()
