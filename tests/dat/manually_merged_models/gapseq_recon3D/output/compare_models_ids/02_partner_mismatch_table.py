#!/usr/bin/env python3
"""Compare microbial partner identities for matched host↔microbiome pools.

This script is meant to complement the match/no-match counts by checking *which* microbial
metabolite each matched host pool is connected to.

Definitions (mirrors the notebook logic):
- Manual merge ("ground truth" pool definition): a pool is considered matched if a feeding-compartment
  metabolite (compartment "f") is fed by BOTH a host connector (H_IEX_*) and a microbial connector
  (M_IEX_*). If an IEX reaction contains both a non-proton and a proton in "f", the non-proton pool is
  used. Proton-only IEX reactions define the proton ("h") pool.
- MeMoMe merge: a pool is considered matched if its translation-compartment metabolite ("t") is
  connected to BOTH model1_* and model2_* via TR_* reactions.

Microbial partner identity:
- Manual: for each matched pool, the partner is the canonical ID of the microbial external metabolite
  (M_* in compartment "e") participating in the corresponding M_IEX_* reaction(s).
- MeMoMe: for each matched pool, the partner is the canonical ID of the microbial external metabolite
  (model2_* in compartment "e") participating in the corresponding TR_model2_* reaction(s).

Outputs:
- tp_partner_mapping.csv: all pools matched in BOTH models, with microbial partner(s) on each side and
  a mismatch flag.
- tp_partner_mismatches.csv: subset of rows where microbial partner sets differ.
"""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import cobra


BASE = Path(__file__).resolve().parent
AUTO_MODEL = BASE.parent / "automatically_merged_metamodel.xml"
MANUAL_MODEL = BASE.parent / "merged_model_2025_prefixed_normalized_diet_fixIEX.xml"


def _load_canonicalizers():
    src = BASE / "00_compare_models_ids.py"
    spec = importlib.util.spec_from_file_location("cmp", src)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.canonical_met_id


def _join(values: set[str]) -> str:
    return "|".join(sorted(values))


def _feed_base_from_manual_iex(rxn, canon_met) -> str | None:
    """Determine which pool base-ID an IEX reaction refers to (feeding compartment 'f').

    Rule: if any non-proton metabolite is present in compartment 'f', use that; otherwise, if the
    reaction is proton-only in 'f', return 'h'.
    """
    feed_f = [m for m in rxn.metabolites if getattr(m, "compartment", None) == "f"]
    if not feed_f:
        return None

    feed_bases = [(m, canon_met(m.id)) for m in feed_f]
    non_proton = [b for _, b in feed_bases if b is not None and b != "h"]
    if non_proton:
        return non_proton[0]

    proton = [b for _, b in feed_bases if b == "h"]
    return "h" if proton else None


def _matched_pools_manual(model, canon_met) -> set[str]:
    host = set()
    micro = set()
    for rxn in model.reactions:
        if rxn.id.startswith("H_IEX_"):
            b = _feed_base_from_manual_iex(rxn, canon_met)
            if b:
                host.add(b)
        elif rxn.id.startswith("M_IEX_"):
            b = _feed_base_from_manual_iex(rxn, canon_met)
            if b:
                micro.add(b)
    return host & micro


def _matched_pools_auto(model, canon_met) -> set[str]:
    matched: set[str] = set()
    for tm in model.metabolites:
        if getattr(tm, "compartment", None) != "t":
            continue
        b = canon_met(tm.id)
        if not b:
            continue
        sides: set[str] = set()
        for rxn in tm.reactions:
            if not rxn.id.startswith("TR_"):
                continue
            for other in rxn.metabolites:
                if other.id == tm.id:
                    continue
                if other.id.startswith("model1_"):
                    sides.add("model1")
                elif other.id.startswith("model2_"):
                    sides.add("model2")
        if sides == {"model1", "model2"}:
            matched.add(b)
    return matched


def _micro_partner_map_manual(model, canon_met) -> dict[str, dict[str, set[str]]]:
    """Map pool base-ID -> partner info from M_IEX_* reactions."""
    out: dict[str, dict[str, set[str]]] = {}
    for rxn in model.reactions:
        if not rxn.id.startswith("M_IEX_"):
            continue
        pool = _feed_base_from_manual_iex(rxn, canon_met)
        if pool is None:
            continue

        info = out.setdefault(
            pool,
            {
                "partner_base_ids": set(),
                "partner_met_ids": set(),
                "partner_names": set(),
                "link_rxn_ids": set(),
            },
        )
        info["link_rxn_ids"].add(rxn.id)
        for met in rxn.metabolites:
            if getattr(met, "compartment", None) == "e" and met.id.startswith("M_"):
                info["partner_base_ids"].add(canon_met(met.id))
                info["partner_met_ids"].add(met.id)
                if met.name:
                    info["partner_names"].add(met.name)
    return out


def _micro_partner_map_auto(model, canon_met) -> dict[str, dict[str, set[str]]]:
    """Map pool base-ID -> partner info from TR_model2_* reactions."""
    out: dict[str, dict[str, set[str]]] = {}
    for rxn in model.reactions:
        if not rxn.id.startswith("TR_model2_"):
            continue

        # This TR reaction may connect one model2 metabolite to one or more translation metabolites.
        pools = {
            canon_met(m.id)
            for m in rxn.metabolites
            if getattr(m, "compartment", None) == "t"
        }
        if not pools:
            continue

        partners = [
            m
            for m in rxn.metabolites
            if getattr(m, "compartment", None) == "e" and m.id.startswith("model2_")
        ]

        for pool in pools:
            if not pool:
                continue
            info = out.setdefault(
                pool,
                {
                    "partner_base_ids": set(),
                    "partner_met_ids": set(),
                    "partner_names": set(),
                    "link_rxn_ids": set(),
                },
            )
            info["link_rxn_ids"].add(rxn.id)
            for met in partners:
                info["partner_base_ids"].add(canon_met(met.id))
                info["partner_met_ids"].add(met.id)
                if met.name:
                    info["partner_names"].add(met.name)
    return out


def main() -> None:
    canon_met = _load_canonicalizers()

    auto = cobra.io.read_sbml_model(str(AUTO_MODEL))
    manual = cobra.io.read_sbml_model(str(MANUAL_MODEL))

    matched_manual = _matched_pools_manual(manual, canon_met)
    matched_auto = _matched_pools_auto(auto, canon_met)
    tp = sorted(matched_manual & matched_auto)

    manual_partner = _micro_partner_map_manual(manual, canon_met)
    auto_partner = _micro_partner_map_auto(auto, canon_met)

    mapping_path = BASE / "tp_partner_mapping.csv"
    mismatch_path = BASE / "tp_partner_mismatches.csv"

    rows = []
    for pool in tp:
        man = manual_partner.get(pool, {})
        au = auto_partner.get(pool, {})
        man_set = set(man.get("partner_base_ids", set()))
        au_set = set(au.get("partner_base_ids", set()))
        mismatch = man_set != au_set
        rows.append(
            {
                "pool_base_id": pool,
                "manual_micro_partner_base_ids": _join(man_set),
                "auto_micro_partner_base_ids": _join(au_set),
                "partner_mismatch": int(mismatch),
                "manual_M_IEX_rxn_ids": _join(set(man.get("link_rxn_ids", set()))),
                "auto_TR_model2_rxn_ids": _join(set(au.get("link_rxn_ids", set()))),
                "manual_micro_met_ids": _join(set(man.get("partner_met_ids", set()))),
                "auto_micro_met_ids": _join(set(au.get("partner_met_ids", set()))),
                "manual_micro_met_names": _join(set(man.get("partner_names", set()))),
                "auto_micro_met_names": _join(set(au.get("partner_names", set()))),
            }
        )

    fieldnames = list(rows[0].keys()) if rows else []
    with mapping_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    mismatches = [r for r in rows if r["partner_mismatch"]]
    with mismatch_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(mismatches)

    print(f"manual matched pools: {len(matched_manual)}")
    print(f"auto matched pools: {len(matched_auto)}")
    print(f"TP pools (matched in both): {len(tp)}")
    print(f"partner mismatches among TP: {len(mismatches)}")
    print(f"wrote: {mapping_path.name} | {mismatch_path.name}")


if __name__ == "__main__":
    main()

