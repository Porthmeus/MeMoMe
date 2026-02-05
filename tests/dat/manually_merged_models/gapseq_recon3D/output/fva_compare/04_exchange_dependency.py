#!/usr/bin/env python3
"""
Compute microbiome-dependency on *exchange reactions* (Auto vs Manual).

This is the EX-reaction analogue of `fva_compare/01_run_fva_analysis.py` +
`fva_compare/02_postprocess_fva_analysis.py`, but focused on the diet exchanges in:

- Auto model: translation compartment "t" (EX_* boundary rxns over *_t metabolites)
- Manual model: feeding compartment "f" (EX_* boundary rxns over *[f] metabolites)

We run FVA twice (baseline/open microbiome vs closed microbiome), and derive:
- dependency_range = (range_open - range_closed) / range_open
- dependency_min   = (min_open   - min_closed)   / min_open
- dependency_max   = (max_open   - max_closed)   / max_open

All comparisons are keyed by *canonical metabolite id* (seed/VMH-normalized),
so auto/manual EX IDs can differ but still compare if they exchange the same metabolite.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis

from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix


# -----------------------------------------------------------------------------
# Inputs / outputs
# -----------------------------------------------------------------------------
ROOT = Path("tests/dat/manually_merged_models/gapseq_recon3D/output")
AUTO_MODEL_PATH = ROOT / "automatically_merged_metamodel.xml"
MANUAL_MODEL_PATH = ROOT / "merged_model_2025_prefixed_normalized_diet_fixIEX.xml"


# -----------------------------------------------------------------------------
# Constraints (kept consistent with host-reaction analysis in 01_run_fva_analysis.py)
# -----------------------------------------------------------------------------
HOST_BIOMASS_LOWER_BOUND = 0.0001
BACTERIAL_BIOMASS_LOWER_BOUND = 0.001

# Force microbiome internal O2 uptake lower bound to 0 (baseline + closed)
AUTO_FORCE_ZERO_O2_RXN = "TR_model2_cpd00007[e]"
MANUAL_FORCE_ZERO_O2_RXN = "M_IEX_cpd00007[e]"


# -----------------------------------------------------------------------------
# Model-specific IDs / prefixes (same idea as 01_run_fva_analysis.py)
# -----------------------------------------------------------------------------
AUTO = {
    "label": "auto",
    "model_path": AUTO_MODEL_PATH,
    "feeding_compartment": "t",
    "host_biomass_id": "model1_biomass_reaction",
    "bacterial_biomass_ids": [
        "model2_Bacterial_Gram_negative_biomass_reaction",
        "model2_Bacterial_Gram_positive_biomass_reaction",
    ],
    "host_connector_prefix": "TR_model1_",
    "microbe_connector_prefix": "TR_model2_",
    "force_zero_o2_rxn": AUTO_FORCE_ZERO_O2_RXN,
}

MANUAL = {
    "label": "manual",
    "model_path": MANUAL_MODEL_PATH,
    "feeding_compartment": "f",
    "host_biomass_id": "H_biomass_reaction",
    "bacterial_biomass_ids": [
        "M_Bacterial_Gram_negative_biomass_reaction",
        "M_Bacterial_Gram_positive_biomass_reaction",
    ],
    "host_connector_prefix": "H_IEX_",
    "microbe_connector_prefix": "M_IEX_",
    "force_zero_o2_rxn": MANUAL_FORCE_ZERO_O2_RXN,
}


# -----------------------------------------------------------------------------
# ID canonicalization (reused concept from compare_models_ids/00_compare_models_ids.py)
# -----------------------------------------------------------------------------
def _strip_model_prefix(value: str) -> str:
    for prefix in ("model1_", "model2_", "H_", "M_"):
        if value.startswith(prefix):
            return value[len(prefix) :]
    return value


def canonical_met_id(met_id: str) -> str:
    base = _strip_model_prefix(met_id)
    cleaned = handle_metabolites_prefix_suffix(base)
    return cleaned if cleaned is not None else base


# -----------------------------------------------------------------------------
# Exchange reaction selection (robust: use model.boundary, not model.exchanges)
# -----------------------------------------------------------------------------
def exchange_boundary_rxns(model: cobra.Model, *, feeding_compartment: str) -> list[cobra.Reaction]:
    """EX_* boundary reactions over a single metabolite in the feeding compartment."""
    rxns: list[cobra.Reaction] = []
    for rxn in model.boundary:
        if not rxn.id.startswith("EX_"):
            continue
        if len(rxn.metabolites) != 1:
            continue
        met = next(iter(rxn.metabolites))
        if met.compartment != feeding_compartment:
            continue
        rxns.append(rxn)
    return rxns


# -----------------------------------------------------------------------------
# Constraint helpers (copied from 01_run_fva_analysis.py semantics)
# -----------------------------------------------------------------------------
def apply_biomass_constraints(model: cobra.Model, *, host_biomass_id: str, bacterial_biomass_ids: list[str], bacterial_lb: float) -> None:
    model.reactions.get_by_id(host_biomass_id).lower_bound = HOST_BIOMASS_LOWER_BOUND
    for rid in bacterial_biomass_ids:
        model.reactions.get_by_id(rid).lower_bound = bacterial_lb


def force_microbiome_o2_zero(model: cobra.Model, rxn_id: str) -> None:
    try:
        model.reactions.get_by_id(rxn_id).lower_bound = 0.0
    except KeyError:
        # Non-fatal: model variants may not contain this reaction.
        pass


def close_microbiome_connectors(model: cobra.Model, *, microbe_connector_prefix: str) -> int:
    """Set bounds to 0 for microbial connector reactions (TR_model2_* or M_IEX_*)."""
    closed = 0
    for rxn in model.reactions:
        if rxn.id.startswith(microbe_connector_prefix):
            rxn.lower_bound = 0.0
            rxn.upper_bound = 0.0
            closed += 1
    return closed


# -----------------------------------------------------------------------------
# FVA + dependency computation (same pattern as 01_run_fva_analysis.py)
# -----------------------------------------------------------------------------
def run_fva(model: cobra.Model, rxn_ids: list[str], *, processes: int | None) -> pd.DataFrame:
    df = flux_variability_analysis(
        model,
        reaction_list=rxn_ids,
        fraction_of_optimum=0.0,
        processes=processes,
    )
    df = df.round(6)
    df["range"] = df["maximum"] - df["minimum"]
    return df


def compute_dependency_table(baseline: pd.DataFrame, closed: pd.DataFrame) -> pd.DataFrame:
    """Return baseline/closed stats + dependency_range/min/max, indexed by rxn id."""
    common = baseline.index.intersection(closed.index)
    base = baseline.loc[common]
    clo = closed.loc[common]

    out = pd.DataFrame(
        {
            "baseline_min": base["minimum"].astype(float),
            "baseline_max": base["maximum"].astype(float),
            "baseline_range": base["range"].astype(float),
            "closed_min": clo["minimum"].astype(float),
            "closed_max": clo["maximum"].astype(float),
            "closed_range": clo["range"].astype(float),
        },
        index=common.astype(str),
    )

    def _rel_change(a: pd.Series, b: pd.Series) -> pd.Series:
        denom = a.to_numpy(dtype=float)
        numer = (a - b).to_numpy(dtype=float)
        res = np.full_like(denom, np.nan, dtype=float)
        mask = denom != 0.0
        res[mask] = numer[mask] / denom[mask]
        return pd.Series(res, index=a.index)

    out["dependency"] = _rel_change(out["baseline_range"], out["closed_range"])
    out["dependency_min"] = _rel_change(out["baseline_min"], out["closed_min"])
    out["dependency_max"] = _rel_change(out["baseline_max"], out["closed_max"])
    return out


def build_exchange_table(
    model: cobra.Model,
    *,
    feeding_compartment: str,
    baseline_fva: pd.DataFrame,
    closed_fva: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build per-exchange table keyed by canonical metabolite id.

    Columns include:
    - exchange_rxn_id, met_id, met_name, exchange_lb/ub
    - baseline/closed min/max/range from FVA
    - dependency_range/min/max
    """
    ex_rxns = exchange_boundary_rxns(model, feeding_compartment=feeding_compartment)

    meta_rows = []
    for rxn in ex_rxns:
        met = next(iter(rxn.metabolites))
        meta_rows.append(
            {
                "exchange_rxn_id": rxn.id,
                "canonical_met_id": canonical_met_id(met.id),
                "met_id": met.id,
                "met_name": met.name,
                "exchange_lb": rxn.lower_bound,
                "exchange_ub": rxn.upper_bound,
            }
        )
    meta = pd.DataFrame(meta_rows).set_index("exchange_rxn_id", drop=True)

    stats = compute_dependency_table(baseline_fva, closed_fva)
    stats.index.name = "exchange_rxn_id"

    df = meta.join(stats, how="inner")
    # Collapse to canonical_met_id (exchange IDs differ between models).
    dup = df["canonical_met_id"].duplicated(keep=False)
    dup_count = int(dup.sum())
    if dup_count:
        print(f"WARNING: {dup_count} duplicated canonical_met_id rows; keeping the first per id.")
    df = df.sort_values(["canonical_met_id", "exchange_rxn_id"]).drop_duplicates(subset=["canonical_met_id"], keep="first")
    df = df.set_index("canonical_met_id", drop=True)
    return df


def build_exchange_connectivity_table(model: cobra.Model, *, feeding_compartment: str, host_connector_prefix: str, microbe_connector_prefix: str) -> pd.DataFrame:
    """For each exchanged metabolite, record whether host/microbe connectors exist."""
    ex_rxns = exchange_boundary_rxns(model, feeding_compartment=feeding_compartment)

    info: dict[str, dict] = {}
    for rxn in ex_rxns:
        met = next(iter(rxn.metabolites))
        cid = canonical_met_id(met.id)
        info[cid] = {
            "canonical_met_id": cid,
            "met_id": met.id,
            "met_name": met.name,
            "exchange_rxn_id": rxn.id,
            "host_connector_rxns": [],
            "microbe_connector_rxns": [],
        }

    def _attach(rxn: cobra.Reaction, key: str) -> None:
        feeding = [m for m in rxn.metabolites if m.compartment == feeding_compartment]
        if len(feeding) != 1:
            return
        cid = canonical_met_id(feeding[0].id)
        info.setdefault(
            cid,
            {
                "canonical_met_id": cid,
                "met_id": None,
                "met_name": None,
                "exchange_rxn_id": None,
                "host_connector_rxns": [],
                "microbe_connector_rxns": [],
            },
        )
        info[cid][key].append(rxn.id)

    for rxn in model.reactions:
        if rxn.id.startswith(host_connector_prefix):
            _attach(rxn, "host_connector_rxns")
        if rxn.id.startswith(microbe_connector_prefix):
            _attach(rxn, "microbe_connector_rxns")

    rows = []
    for cid, d in info.items():
        host_rxns = sorted(d["host_connector_rxns"])
        micro_rxns = sorted(d["microbe_connector_rxns"])
        rows.append(
            {
                "canonical_met_id": cid,
                "met_id": d.get("met_id"),
                "met_name": d.get("met_name"),
                "exchange_rxn_id": d.get("exchange_rxn_id"),
                "has_host_connector": bool(host_rxns),
                "has_microbe_connector": bool(micro_rxns),
                "host_connector_count": len(host_rxns),
                "microbe_connector_count": len(micro_rxns),
                "host_connector_rxns": "|".join(host_rxns),
                "microbe_connector_rxns": "|".join(micro_rxns),
            }
        )
    return pd.DataFrame(rows).sort_values("canonical_met_id")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "fva_compare" / "exchange_analysis",
        help="Directory to write exchange-FVA and dependency outputs.",
    )
    parser.add_argument("--processes", type=int, default=None, help="FVA multiprocessing workers.")
    parser.add_argument("--force", action="store_true", help="Recompute FVA even if cached CSVs exist.")
    args = parser.parse_args()

    outdir: Path = args.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    auto_model = cobra.io.read_sbml_model(str(AUTO["model_path"]))
    manual_model = cobra.io.read_sbml_model(str(MANUAL["model_path"]))

    auto_ex_ids = [r.id for r in exchange_boundary_rxns(auto_model, feeding_compartment=AUTO["feeding_compartment"])]
    manual_ex_ids = [r.id for r in exchange_boundary_rxns(manual_model, feeding_compartment=MANUAL["feeding_compartment"])]
    print(f"auto EX in '{AUTO['feeding_compartment']}': {len(auto_ex_ids)}")
    print(f"manual EX in '{MANUAL['feeding_compartment']}': {len(manual_ex_ids)}")

    # ----- Baseline (microbiome open)
    apply_biomass_constraints(
        auto_model,
        host_biomass_id=AUTO["host_biomass_id"],
        bacterial_biomass_ids=AUTO["bacterial_biomass_ids"],
        bacterial_lb=BACTERIAL_BIOMASS_LOWER_BOUND,
    )
    apply_biomass_constraints(
        manual_model,
        host_biomass_id=MANUAL["host_biomass_id"],
        bacterial_biomass_ids=MANUAL["bacterial_biomass_ids"],
        bacterial_lb=BACTERIAL_BIOMASS_LOWER_BOUND,
    )
    force_microbiome_o2_zero(auto_model, AUTO["force_zero_o2_rxn"])
    force_microbiome_o2_zero(manual_model, MANUAL["force_zero_o2_rxn"])

    auto_base_path = outdir / "auto_exchange_fva_baseline.csv"
    manual_base_path = outdir / "manual_exchange_fva_baseline.csv"
    if auto_base_path.exists() and manual_base_path.exists() and not args.force:
        auto_base = pd.read_csv(auto_base_path, index_col=0)
        manual_base = pd.read_csv(manual_base_path, index_col=0)
    else:
        auto_base = run_fva(auto_model, auto_ex_ids, processes=args.processes)
        manual_base = run_fva(manual_model, manual_ex_ids, processes=args.processes)
        auto_base.to_csv(auto_base_path)
        manual_base.to_csv(manual_base_path)

    # ----- Closed microbiome (microbial connectors forced to 0)
    auto_closed = auto_model.copy()
    manual_closed = manual_model.copy()
    closed_auto_n = close_microbiome_connectors(auto_closed, microbe_connector_prefix=AUTO["microbe_connector_prefix"])
    closed_manual_n = close_microbiome_connectors(manual_closed, microbe_connector_prefix=MANUAL["microbe_connector_prefix"])
    print(f"closed connectors: auto={closed_auto_n} manual={closed_manual_n}")

    # Remove bacterial biomass constraints to keep the models feasible after closure.
    apply_biomass_constraints(
        auto_closed,
        host_biomass_id=AUTO["host_biomass_id"],
        bacterial_biomass_ids=AUTO["bacterial_biomass_ids"],
        bacterial_lb=0.0,
    )
    apply_biomass_constraints(
        manual_closed,
        host_biomass_id=MANUAL["host_biomass_id"],
        bacterial_biomass_ids=MANUAL["bacterial_biomass_ids"],
        bacterial_lb=0.0,
    )
    force_microbiome_o2_zero(auto_closed, AUTO["force_zero_o2_rxn"])
    force_microbiome_o2_zero(manual_closed, MANUAL["force_zero_o2_rxn"])

    auto_closed_path = outdir / "auto_exchange_fva_closed.csv"
    manual_closed_path = outdir / "manual_exchange_fva_closed.csv"
    if auto_closed_path.exists() and manual_closed_path.exists() and not args.force:
        auto_clo = pd.read_csv(auto_closed_path, index_col=0)
        manual_clo = pd.read_csv(manual_closed_path, index_col=0)
    else:
        auto_clo = run_fva(auto_closed, auto_ex_ids, processes=args.processes)
        manual_clo = run_fva(manual_closed, manual_ex_ids, processes=args.processes)
        auto_clo.to_csv(auto_closed_path)
        manual_clo.to_csv(manual_closed_path)

    # ----- Build per-model exchange tables (keyed by canonical metabolite id)
    auto_dep = build_exchange_table(
        auto_model,
        feeding_compartment=AUTO["feeding_compartment"],
        baseline_fva=auto_base,
        closed_fva=auto_clo,
    )
    manual_dep = build_exchange_table(
        manual_model,
        feeding_compartment=MANUAL["feeding_compartment"],
        baseline_fva=manual_base,
        closed_fva=manual_clo,
    )

    # Connector presence tables (host vs microbiome)
    auto_conn = build_exchange_connectivity_table(
        auto_model,
        feeding_compartment=AUTO["feeding_compartment"],
        host_connector_prefix=AUTO["host_connector_prefix"],
        microbe_connector_prefix=AUTO["microbe_connector_prefix"],
    )
    manual_conn = build_exchange_connectivity_table(
        manual_model,
        feeding_compartment=MANUAL["feeding_compartment"],
        host_connector_prefix=MANUAL["host_connector_prefix"],
        microbe_connector_prefix=MANUAL["microbe_connector_prefix"],
    )
    auto_conn.to_csv(outdir / "auto_exchange_connectivity.csv", index=False)
    manual_conn.to_csv(outdir / "manual_exchange_connectivity.csv", index=False)

    # Restrict exchange analysis to metabolites connected to BOTH host and microbiome.
    # This drops host-only or microbe-only "diet" metabolites that cannot mediate
    # host <-> microbiome interactions via the feeding/translation compartment.
    auto_both_connected = set(
        auto_conn.loc[auto_conn["has_host_connector"] & auto_conn["has_microbe_connector"], "canonical_met_id"].astype(str)
    )
    manual_both_connected = set(
        manual_conn.loc[manual_conn["has_host_connector"] & manual_conn["has_microbe_connector"], "canonical_met_id"].astype(str)
    )
    auto_both_connected &= set(auto_dep.index.astype(str))
    manual_both_connected &= set(manual_dep.index.astype(str))
    both_connected_overlap = auto_both_connected & manual_both_connected

    print(f"auto both-connected exchanged metabolites: {len(auto_both_connected)}")
    print(f"manual both-connected exchanged metabolites: {len(manual_both_connected)}")
    print(f"both-connected overlap (auto ∩ manual): {len(both_connected_overlap)}")

    auto_dep = auto_dep.loc[sorted(auto_both_connected)]
    manual_dep = manual_dep.loc[sorted(manual_both_connected)]

    auto_dep.to_csv(outdir / "auto_exchange_dependency.csv")
    manual_dep.to_csv(outdir / "manual_exchange_dependency.csv")

    # ----- Intersection-only comparison table (canonical_met_id index)
    common = auto_dep.index.intersection(manual_dep.index)
    comparison = pd.DataFrame(index=common)
    for prefix, src in (("auto", auto_dep), ("manual", manual_dep)):
        comparison[f"baseline_range_{prefix}"] = src.loc[common, "baseline_range"]
        comparison[f"closed_range_{prefix}"] = src.loc[common, "closed_range"]
        comparison[f"dependency_{prefix}"] = src.loc[common, "dependency"]
        comparison[f"baseline_min_{prefix}"] = src.loc[common, "baseline_min"]
        comparison[f"closed_min_{prefix}"] = src.loc[common, "closed_min"]
        comparison[f"dependency_min_{prefix}"] = src.loc[common, "dependency_min"]
        comparison[f"baseline_max_{prefix}"] = src.loc[common, "baseline_max"]
        comparison[f"closed_max_{prefix}"] = src.loc[common, "closed_max"]
        comparison[f"dependency_max_{prefix}"] = src.loc[common, "dependency_max"]

    comparison["delta_dependency"] = comparison["dependency_auto"] - comparison["dependency_manual"]
    comparison["abs_delta_dependency"] = comparison["delta_dependency"].abs()
    comparison["delta_dependency_min"] = comparison["dependency_min_auto"] - comparison["dependency_min_manual"]
    comparison["abs_delta_dependency_min"] = comparison["delta_dependency_min"].abs()
    comparison["delta_dependency_max"] = comparison["dependency_max_auto"] - comparison["dependency_max_manual"]
    comparison["abs_delta_dependency_max"] = comparison["delta_dependency_max"].abs()
    comparison.to_csv(outdir / "auto_manual_exchange_dependency.csv")

    # ----- Convenience outputs (top dependent exchanges + top deltas)
    results_dir = outdir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    auto_dep.sort_values("dependency", ascending=False).head(200).to_csv(results_dir / "auto_top_dependent_exchanges.csv")
    manual_dep.sort_values("dependency", ascending=False).head(200).to_csv(results_dir / "manual_top_dependent_exchanges.csv")
    comparison.sort_values("abs_delta_dependency", ascending=False).head(200).to_csv(results_dir / "top_abs_delta_exchange_dependency.csv")

    auto_dep.dropna(subset=["dependency_min"]).sort_values("dependency_min", ascending=False).head(200).to_csv(
        results_dir / "auto_top_dependent_exchanges_by_min.csv"
    )
    manual_dep.dropna(subset=["dependency_min"]).sort_values("dependency_min", ascending=False).head(200).to_csv(
        results_dir / "manual_top_dependent_exchanges_by_min.csv"
    )
    comparison.dropna(subset=["abs_delta_dependency_min"]).sort_values("abs_delta_dependency_min", ascending=False).head(200).to_csv(
        results_dir / "top_abs_delta_exchange_dependency_min.csv"
    )

    auto_dep.dropna(subset=["dependency_max"]).sort_values("dependency_max", ascending=False).head(200).to_csv(
        results_dir / "auto_top_dependent_exchanges_by_max.csv"
    )
    manual_dep.dropna(subset=["dependency_max"]).sort_values("dependency_max", ascending=False).head(200).to_csv(
        results_dir / "manual_top_dependent_exchanges_by_max.csv"
    )
    comparison.dropna(subset=["abs_delta_dependency_max"]).sort_values("abs_delta_dependency_max", ascending=False).head(200).to_csv(
        results_dir / "top_abs_delta_exchange_dependency_max.csv"
    )

    print(f"exchange dependency intersection size: {len(common)}")
    print(f"outputs written to: {outdir}")


if __name__ == "__main__":
    main()
