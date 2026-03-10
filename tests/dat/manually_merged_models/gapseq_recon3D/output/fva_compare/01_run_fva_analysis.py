#!/usr/bin/env python3
"""Run FVAs and save raw results (no postprocessing).
Intended to compare FVA results between automatically merged and manually merged models.
"""

# 1) Imports
from pathlib import Path

import cobra
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis


ROOT = Path("tests/dat/manually_merged_models/gapseq_recon3D/output")
AUTO_MODEL_PATH = ROOT / "automatically_merged_metamodel.xml"
MANUAL_MODEL_PATH = ROOT / "merged_model_2025_prefixed_normalized_diet_fixIEX.xml"
OUT_DIR = ROOT / "fva_compare"
OUT_DIR.mkdir(exist_ok=True, parents=True)

# Biomass reaction IDs
AUTO_HOST_BIOMASS_ID = "model1_biomass_reaction"
AUTO_BACTERIAL_BIOMASS_IDS = ["model2_Bacterial_Gram_negative_biomass_reaction", "model2_Bacterial_Gram_positive_biomass_reaction"]
MANUAL_HOST_BIOMASS_ID = "H_biomass_reaction"
MANUAL_BACTERIAL_BIOMASS_IDS = ["M_Bacterial_Gram_negative_biomass_reaction", "M_Bacterial_Gram_positive_biomass_reaction"]
#BIOMASS LOWER BOUNDS
HOST_BIOMASS_LOWER_BOUND = 0.0001
BACTERIAL_BIOMASS_LOWER_BOUND = 0.001
#Create a dictionaries for easy access to biomass IDs and model prefixes
BIOMASS_IDS = {
    "auto": {
        "host": [AUTO_HOST_BIOMASS_ID],
        "bacterial": AUTO_BACTERIAL_BIOMASS_IDS
    },
    "manual": {
        "host": [MANUAL_HOST_BIOMASS_ID],
        "bacterial": MANUAL_BACTERIAL_BIOMASS_IDS
    }
}

SUBMODEL_PREFIXES = {
    "auto": {"host":"model1_", "bacterial":"model2_"},
    "manual": {"host":"H_", "bacterial":"M_"}
}

#Load models
auto_model = cobra.io.read_sbml_model(str(AUTO_MODEL_PATH))
manual_model = cobra.io.read_sbml_model(str(MANUAL_MODEL_PATH))

#Choose host reaction list
AUTO_HOST_RXN_PREFIXES = (
    "model1_",
    # demand reactions in the auto model are encoded as DM_model1_*
    "DM_model1_",
    # sink reactions in the auto model are encoded as sink_model1_*
    "sink_model1_",
)

def is_auto_host_reaction_id(rxn_id: str) -> bool:
    return rxn_id.startswith(AUTO_HOST_RXN_PREFIXES)

def is_manual_host_reaction_id(rxn_id: str) -> bool:
    # Manual host reactions are prefixed with H_ (including H_DM_* and H_sink_*).
    return rxn_id.startswith("H_")

auto_host_rxns = [rxn.id for rxn in auto_model.reactions if is_auto_host_reaction_id(rxn.id)]
manual_host_rxns = [rxn.id for rxn in manual_model.reactions if is_manual_host_reaction_id(rxn.id)]

#Apply host-only constraints (same for auto + manual)
#   - set HOST_BIOMASS_ID lower bound to HOST_BIOMASS_LOWER_BOUND
auto_model.reactions.get_by_id(AUTO_HOST_BIOMASS_ID).lower_bound = HOST_BIOMASS_LOWER_BOUND
manual_model.reactions.get_by_id(MANUAL_HOST_BIOMASS_ID).lower_bound = HOST_BIOMASS_LOWER_BOUND
#   - set each bacterial biomass lower bound to BACTERIAL_BIOMASS_LOWER_BOUND
for b_biomass_id in AUTO_BACTERIAL_BIOMASS_IDS:
    auto_model.reactions.get_by_id(b_biomass_id).lower_bound = BACTERIAL_BIOMASS_LOWER_BOUND
for b_biomass_id in MANUAL_BACTERIAL_BIOMASS_IDS:
    manual_model.reactions.get_by_id(b_biomass_id).lower_bound = BACTERIAL_BIOMASS_LOWER_BOUND

#Force internal oxygen uptake of the microbiome to 0
auto_model.reactions.get_by_id("TR_model2_cpd00007[e]").lower_bound = 0.0
manual_model.reactions.get_by_id("M_IEX_cpd00007[e]").lower_bound = 0.0

#Run baseline FVA (microbiome open)
auto_base = flux_variability_analysis(auto_model, reaction_list=auto_host_rxns, fraction_of_optimum=0.0)
manual_base = flux_variability_analysis(manual_model, reaction_list=manual_host_rxns, fraction_of_optimum=0.0)
#round the results to the 6th decimal to avoid numerical instability issues in comparisons
auto_base = auto_base.round(6)
manual_base = manual_base.round(6)
#calculate the range column
auto_base["range"] = auto_base["maximum"] - auto_base["minimum"]
manual_base["range"] = manual_base["maximum"] - manual_base["minimum"]
#Save results
auto_base.to_csv(OUT_DIR / "auto_host_fva_baseline.csv")
manual_base.to_csv(OUT_DIR / "manual_host_fva_baseline.csv")
#


def close_auto_microbiome_exchanges(model):
    #close bacterial exchanges for microbiome in auto model
    for rxn in model.reactions:
        if rxn.id.startswith("TR_model2_"):
            rxn.lower_bound = 0.0
            rxn.upper_bound = 0.0
    # remove bacterial biomass constraints to keep the model feasible
    for b_biomass_id in AUTO_BACTERIAL_BIOMASS_IDS:
        model.reactions.get_by_id(b_biomass_id).lower_bound = 0.0
    

def close_manual_microbiome_exchanges(model):
    #close bacterial exchanges for microbiome in manual model
    for rxn in model.reactions:
        if rxn.id.startswith("M_IEX_"):
            rxn.lower_bound = 0.0
            rxn.upper_bound = 0.0
    # remove bacterial biomass constraints to keep the model feasible
    for b_biomass_id in MANUAL_BACTERIAL_BIOMASS_IDS:
        model.reactions.get_by_id(b_biomass_id).lower_bound = 0.0


# Close microbiome exchange reactions, run fva and save results
auto_closed = auto_model.copy()
manual_closed = manual_model.copy()
close_auto_microbiome_exchanges(auto_closed)
close_manual_microbiome_exchanges(manual_closed)
closed_auto_fva = flux_variability_analysis(auto_closed, reaction_list=auto_host_rxns, fraction_of_optimum=0.0)
closed_manual_fva = flux_variability_analysis(manual_closed, reaction_list=manual_host_rxns, fraction_of_optimum=0.0)
#round the results to the 6th decimal to avoid numerical instability issues in comparisons
closed_auto_fva = closed_auto_fva.round(6)
closed_manual_fva = closed_manual_fva.round(6)
#calculate the range column
closed_auto_fva["range"] = closed_auto_fva["maximum"] - closed_auto_fva["minimum"]
closed_manual_fva["range"] = closed_manual_fva["maximum"] - closed_manual_fva["minimum"]
#save results
closed_auto_fva.to_csv(OUT_DIR / "auto_host_fva_closed.csv")
closed_manual_fva.to_csv(OUT_DIR / "manual_host_fva_closed.csv")

