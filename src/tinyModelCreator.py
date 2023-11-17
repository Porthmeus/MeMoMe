import cobra
import pandas as pd
import os

#HARDCODED PARAMETERS
model_file_path = "../tests/dat/e_coli_vmh.xml"
biomass_id = "biomass525" # biomass id of the model

desired_biomass_components =['biomass[c]',
                             'h[c]',
                             'pi[c]',
                             'ppi[c]',
                             'ala_L[c]',
                             'amet[c]',
                             'arg_L[c]',
                             'asn_L[c]',
                             'asp_L[c]',
                             'atp[c]',
                             'ca2[c]',
                             'cl[c]',
                             'proteinsynth[c]',
                             'dnarep[c]',
                             'rnatrans[c]'] #arbitrarily chosen from the biomass components

required_excretions = ['EX_ala_D(e)','EX_ala_L(e)','EX_trp_L(e)',
 'EX_co2(e)','EX_h2o(e)',
 'EX_etoh(e)', 'EX_ac(e)','EX_h2(e)','EX_for(e)','EX_h(e)',
 'EX_lac_D(e)','EX_lac_L(e)'] #manually chosen from the list of potentially excretable compounds (fva results)

desired_medium = dict() #setting the mandatory components of the final medium and their constraints
desired_medium["EX_adn(e)"] = 1.0
desired_medium["EX_gsn(e)"] = 1.0
desired_medium["EX_nmn(e)"] = 100.0
desired_medium["EX_glc_D(e)"] = 1.0
desired_medium["EX_h(e)"] = 10.0
desired_medium["EX_pi(e)"] = 10.0
desired_medium["EX_h2o(e)"] = 1000.0
desired_medium["EX_o2(e)"] = 10.0

new_model_id = "mock_e_coli_vmh"
new_model_name = "Simplified model obtained subsetting vmh's ecoli K12_substr_MG1655"
output_file = "../tests/dat/subset_of_e_coli_vmh.xml"

#MAIN SCRIPT
#load the full e_coli model
model = cobra.io.read_sbml_model(model_file_path)
model

#keep copy of the old medium and old biomass function
old_medium = model.medium.copy()
biomass_func = model.reactions.get_by_id(biomass_id)
old_biomass = biomass_func.copy()

#create a reduced version of the biomass reaction
reactions_toadd = cobra.core.DictList()
new_react = cobra.core.reaction.Reaction("biomass_simplified")
new_react.name = "biomass with less metabolic requirement"
new_react.lower_bound = 0
new_react.upper_bound = 1000
new_react.reversibility = False
new_react_mets = dict()
for met in desired_biomass_components:
    new_react_met = model.metabolites.get_by_id(met)
    new_react_mets[new_react_met] = model.reactions.get_by_id(biomass_id).metabolites[model.metabolites.get_by_id(met)]
new_react.add_metabolites(new_react_mets)
reactions_toadd.add(new_react)
model.add_reactions(reactions_toadd)

#remove old biomass reaction
model.remove_reactions(cobra.core.DictList([model.reactions.get_by_id(biomass_id)]))

#set the new biomass reaction as objective
model.objective = "biomass_simplified"

#calculate optimal growth with new biomass function
max_growth = model.slim_optimize()

#obtain a minimal medium able to support at least 10% of the maximal biomass growth rate
required_medium = cobra.medium.minimal_medium(model, max_growth/10)

#close all the medium's influxes
for reac_id in model.medium.keys():
    model.reactions.get_by_id(reac_id).lower_bound = 0.0

#include the medium components specified by the user
for reac_id in desired_medium.keys():
    required_medium[reac_id] = desired_medium[reac_id]

#set the new medium
for reac_id in required_medium.keys():
    model.reactions.get_by_id(reac_id).lower_bound = - required_medium[reac_id]

#use fva to calculate the excretion and uptakes that the model can achieve while keeping 10% optimality
fva_sol = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0.1)
active_fva_ids = list(set(fva_sol[fva_sol["minimum"] < -1e-06].index)
                      .union(set(fva_sol[fva_sol["maximum"] > 1e-06].index)))
active_fva_exchanges = sorted([r_id for r_id in active_fva_ids if "EX_" in r_id])
possible_excretions = sorted([r_id for r_id in active_fva_exchanges if fva_sol.loc[r_id].maximum > 1e-06])
possible_uptakes = sorted([r_id for r_id in active_fva_exchanges if fva_sol.loc[r_id].minimum < -1e-06])

#select all the sinks and demands that can carry a flux
boundary_ids = {r.id for r in model.boundary}
other_reacs_like_sinks_or_demands = sorted([r_id for r_id in boundary_ids.intersection(set(active_fva_ids)) 
                                              if r_id not in set(possible_excretions).union(set(possible_uptakes))])

#remove from the model all the boundary reactions that are not feasible uptakes, excretions, sinks or demands
reacs_to_remove = list()
for r in model.boundary:
    if r.id not in set(possible_excretions).union(set(possible_uptakes)).union(other_reacs_like_sinks_or_demands):
        reacs_to_remove.append(r)        
model.remove_reactions(reacs_to_remove)

#obtain a small list of reactions that support a 10% maximal biomass production(not proven to be minimal, but with 5 iterations, and using pFBA,the total number of reactions should be small). The excretion of all the compounds required by the user will be guaranteed. All medium reactions required by the user will be kept
for i in range(5):
    model.reactions.get_by_id("biomass_simplified").lower_bound = 0.1 * model.slim_optimize()
    reactions_to_keep = set(required_excretions).union(set(required_medium.keys()))
    for reac in required_excretions:
        model.objective = reac
        pfba_sol = cobra.flux_analysis.pfba(model)
        pfba_sol = pfba_sol.to_frame()
        reactions_to_keep = reactions_to_keep.union(set(pfba_sol[abs(pfba_sol["fluxes"])>1e-06].index))
        if(pfba_sol["fluxes"][reac] <= 1e-06):
            raise Exception("warning: required product " + reac + " cannot be produced by the model")
    model.reactions.get_by_id("biomass_simplified").lower_bound = 0.0
    model.objective = "biomass_simplified"
    reacs_to_remove = list()
    for reac in model.reactions:
        if reac.id not in reactions_to_keep:
            reacs_to_remove.append(reac)
    model.remove_reactions(reacs_to_remove)
    #remove all unused metabolites and genes
    model, removed_metabolites = cobra.manipulation.delete.prune_unused_metabolites(model)
    model_reac_ids = [r.id for r in model.reactions]
    genes_to_delete = list()
    for gene in model.genes:
        useful_gene = False
        for r in gene.reactions:
            if r.id in model_reac_ids:
                useful_gene = True
        if not useful_gene:
            genes_to_delete.append(gene)
    cobra.manipulation.remove_genes(model, gene_list = genes_to_delete)

#set new id and name to the model
model.id = new_model_id
model.name = new_model_name

#save the model
for reac in model.reactions:
    reac.lower_bound = float(reac.lower_bound)
    reac.upper_bound = float(reac.upper_bound)
cobra.io.write_sbml_model(model,output_file)