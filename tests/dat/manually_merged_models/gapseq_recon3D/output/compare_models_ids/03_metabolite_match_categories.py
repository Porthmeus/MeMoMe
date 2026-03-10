#!/usr/bin/env python3
"""Categorize host↔microbiome metabolite matches between manual and MeMoMe merges.

This script assigns each boundary "pool" metabolite (canonical base-ID) to a mutually exclusive
category under the assumption that the manual merge is the ground truth:

Pool-level match definition (mirrors the notebook + partner-mismatch logic):
- Manual merge: a pool is matched if a feeding-compartment metabolite ("f") is fed by BOTH a host
  connector (H_IEX_*) and a microbial connector (M_IEX_*). If an IEX contains both a non-proton and
  a proton in "f", the non-proton metabolite defines the pool. Proton-only IEX reactions define the
  proton ("h") pool.
- MeMoMe merge: a pool is matched if its translation-compartment metabolite ("t") is connected to
  BOTH model1_* and model2_* via TR_* reactions.

Categories (mutually exclusive):
- true_match: matched in both models AND the microbial partner identity matches.
- wrong_match: matched in both models BUT the microbial partner identity differs.
- false_negative: matched in manual, not matched in MeMoMe.
- false_positive: matched in MeMoMe, not matched in manual.
- true_negative: not matched in either model.

Universe:
All canonical pool base-IDs that are referenced by at least one connector in either model:
- Manual: any pool referenced by an H_IEX_* or M_IEX_* connector.
- MeMoMe: any pool referenced by a translation-compartment metabolite ("t") that participates in
  at least one TR_* reaction connecting to model1_* or model2_*.
"""

from __future__ import annotations

import csv
import importlib.util
from collections import Counter
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
    return "|".join(sorted(v for v in values if v))


def _feed_base_from_manual_iex(rxn, canon_met) -> str | None:
    feed_f = [m for m in rxn.metabolites if getattr(m, "compartment", None) == "f"]
    if not feed_f:
        return None

    feed_bases = [(m, canon_met(m.id)) for m in feed_f]
    non_proton = [b for _, b in feed_bases if b is not None and b != "h"]
    if non_proton:
        return non_proton[0]

    proton = [b for _, b in feed_bases if b == "h"]
    return "h" if proton else None


def _manual_pool_index(model, canon_met):
    host_pools: dict[str, set[str]] = {}
    micro_pools: dict[str, set[str]] = {}
    pool_met_ids: dict[str, set[str]] = {}
    pool_names: dict[str, set[str]] = {}

    for rxn in model.reactions:
        if not (rxn.id.startswith("H_IEX_") or rxn.id.startswith("M_IEX_")):
            continue
        pool = _feed_base_from_manual_iex(rxn, canon_met)
        if pool is None:
            continue

        # Record which metabolite in "f" defines this pool (prefer the non-proton base).
        chosen = None
        for met in rxn.metabolites:
            if getattr(met, "compartment", None) != "f":
                continue
            if canon_met(met.id) == pool:
                chosen = met
                break
        if chosen is not None:
            pool_met_ids.setdefault(pool, set()).add(chosen.id)
            if chosen.name:
                pool_names.setdefault(pool, set()).add(chosen.name)

        if rxn.id.startswith("H_IEX_"):
            host_pools.setdefault(pool, set()).add(rxn.id)
        else:
            micro_pools.setdefault(pool, set()).add(rxn.id)

    candidate = set(host_pools) | set(micro_pools)
    matched = set(host_pools) & set(micro_pools)
    return {
        "candidate": candidate,
        "matched": matched,
        "H_IEX_rxn_ids": host_pools,
        "M_IEX_rxn_ids": micro_pools,
        "pool_met_ids": pool_met_ids,
        "pool_names": pool_names,
    }


def _auto_pool_index(model, canon_met):
    side_map: dict[str, set[str]] = {}
    pool_met_ids: dict[str, set[str]] = {}
    pool_names: dict[str, set[str]] = {}
    tr_model1: dict[str, set[str]] = {}
    tr_model2: dict[str, set[str]] = {}

    for tm in model.metabolites:
        if getattr(tm, "compartment", None) != "t":
            continue
        pool = canon_met(tm.id)
        if not pool:
            continue

        pool_met_ids.setdefault(pool, set()).add(tm.id)
        if tm.name:
            pool_names.setdefault(pool, set()).add(tm.name)

        for rxn in tm.reactions:
            if not rxn.id.startswith("TR_"):
                continue

            connects: set[str] = set()
            for other in rxn.metabolites:
                if other.id == tm.id:
                    continue
                if other.id.startswith("model1_"):
                    connects.add("model1")
                elif other.id.startswith("model2_"):
                    connects.add("model2")

            if "model1" in connects:
                tr_model1.setdefault(pool, set()).add(rxn.id)
            if "model2" in connects:
                tr_model2.setdefault(pool, set()).add(rxn.id)

            side_map.setdefault(pool, set()).update(connects)

    candidate = set(side_map)
    matched = {p for p, sides in side_map.items() if sides == {"model1", "model2"}}
    return {
        "candidate": candidate,
        "matched": matched,
        "TR_model1_rxn_ids": tr_model1,
        "TR_model2_rxn_ids": tr_model2,
        "pool_met_ids": pool_met_ids,
        "pool_names": pool_names,
    }


def _micro_partner_map_manual(model, canon_met) -> dict[str, dict[str, set[str]]]:
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
            },
        )
        for met in rxn.metabolites:
            if getattr(met, "compartment", None) == "e" and met.id.startswith("M_"):
                info["partner_base_ids"].add(canon_met(met.id))
                info["partner_met_ids"].add(met.id)
                if met.name:
                    info["partner_names"].add(met.name)
    return out


def _micro_partner_map_auto(model, canon_met) -> dict[str, dict[str, set[str]]]:
    out: dict[str, dict[str, set[str]]] = {}
    for rxn in model.reactions:
        if not rxn.id.startswith("TR_model2_"):
            continue

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
                },
            )
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

    manual_idx = _manual_pool_index(manual, canon_met)
    auto_idx = _auto_pool_index(auto, canon_met)

    manual_partner = _micro_partner_map_manual(manual, canon_met)
    auto_partner = _micro_partner_map_auto(auto, canon_met)

    universe = sorted(manual_idx["candidate"] | auto_idx["candidate"])

    out_path = BASE / "metabolite_match_categories.csv"
    counts_path = BASE / "metabolite_match_category_counts.csv"

    rows = []
    for pool in universe:
        manual_matched = pool in manual_idx["matched"]
        auto_matched = pool in auto_idx["matched"]

        man_p = manual_partner.get(pool, {})
        au_p = auto_partner.get(pool, {})
        man_set = set(man_p.get("partner_base_ids", set()))
        au_set = set(au_p.get("partner_base_ids", set()))

        if manual_matched and auto_matched:
            category = "true_match" if man_set == au_set else "wrong_match"
        elif manual_matched and not auto_matched:
            category = "false_negative"
        elif (not manual_matched) and auto_matched:
            category = "false_positive"
        else:
            category = "true_negative"

        rows.append(
            {
                "pool_base_id": pool,
                "category": category,
                "manual_matched": int(manual_matched),
                "auto_matched": int(auto_matched),
                "manual_micro_partner_base_ids": _join(man_set),
                "auto_micro_partner_base_ids": _join(au_set),
                "manual_micro_partner_names": _join(set(man_p.get("partner_names", set()))),
                "auto_micro_partner_names": _join(set(au_p.get("partner_names", set()))),
                "manual_pool_met_ids": _join(set(manual_idx["pool_met_ids"].get(pool, set()))),
                "auto_pool_met_ids": _join(set(auto_idx["pool_met_ids"].get(pool, set()))),
                "manual_pool_names": _join(set(manual_idx["pool_names"].get(pool, set()))),
                "auto_pool_names": _join(set(auto_idx["pool_names"].get(pool, set()))),
                "manual_H_IEX_rxn_ids": _join(set(manual_idx["H_IEX_rxn_ids"].get(pool, set()))),
                "manual_M_IEX_rxn_ids": _join(set(manual_idx["M_IEX_rxn_ids"].get(pool, set()))),
                "auto_TR_model1_rxn_ids": _join(set(auto_idx["TR_model1_rxn_ids"].get(pool, set()))),
                "auto_TR_model2_rxn_ids": _join(set(auto_idx["TR_model2_rxn_ids"].get(pool, set()))),
            }
        )

    with out_path.open("w", newline="") as f:
        fieldnames = list(rows[0].keys()) if rows else []
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    category_counts = Counter(r["category"] for r in rows)
    with counts_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["category", "count"])
        w.writeheader()
        for cat, n in sorted(category_counts.items()):
            w.writerow({"category": cat, "count": n})

    # Print a concise summary for quick inspection.
    print(f"manual candidates: {len(manual_idx['candidate'])} | matched: {len(manual_idx['matched'])}")
    print(f"auto candidates:   {len(auto_idx['candidate'])} | matched: {len(auto_idx['matched'])}")
    print(f"universe pools:    {len(universe)}")
    for cat in ("true_match", "wrong_match", "false_negative", "false_positive", "true_negative"):
        print(f"{cat}: {category_counts.get(cat, 0)}")
    print(f"wrote: {out_path.name} | {counts_path.name}")


if __name__ == "__main__":
    main()

