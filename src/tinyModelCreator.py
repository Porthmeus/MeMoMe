
# Import necessary libraries
import argparse
import cobra
import pandas as pd
import os

# Setup argparse to handle input arguments
parser = argparse.ArgumentParser(description="Modify a metabolic model based on specified parameters.")
parser.add_argument("--model_file_path", type=str, default="../tests/dat/script_accessories/e_coli_vmh.xml", help="Path to the model file.")
parser.add_argument("--biomass_id", type=str, default="biomass525", help="ID of the biomass reaction in the model.")
parser.add_argument("--desired_biomass_components_file", type=str, default="../tests/dat/script_accessories/desired_biomass_components.txt", help="File containing desired biomass components.")
parser.add_argument("--required_excretions_file", type=str, default="../tests/dat/script_accessories/required_excretions.txt", help="File containing required excretions.")
parser.add_argument("--desired_medium_file", type=str, default="../tests/dat/script_accessories/desired_medium.csv", help="CSV file containing desired medium components and constraints.")
parser.add_argument("--new_model_id", type=str ,default="minimal_model", help="New ID for the modified model.")
parser.add_argument("--new_model_name", type=str, default="minimal_model", help="New name for the modified model.")
parser.add_argument("--output_file", type=str, default="../tests/dat/script_accessories/modified_model.xml", help="Path for the output modified model.")

args = parser.parse_args()

# Read in the desired biomass components
with open(args.desired_biomass_components_file) as file:
    desired_biomass_components = [line.strip() for line in file]

# Read in the required excretions
with open(args.required_excretions_file) as file:
    required_excretions = [line.strip() for line in file]

# Load the desired medium from CSV
desired_medium_df = pd.read_csv(args.desired_medium_file)
desired_medium = pd.Series(desired_medium_df['Constraint'].values, index=desired_medium_df['Component']).to_dict()

# Load the full model
model = cobra.io.read_sbml_model(args.model_file_path)

# Keep a copy of the old medium and old biomass function
old_medium = model.medium.copy()
biomass_func = model.reactions.get_by_id(args.biomass_id)
old_biomass = biomass_func.copy()

# Create a reduced version of the biomass reaction 
reactions_toadd = cobra.core.DictList()
new_react = cobra.core.reaction.Reaction("biomass_simplified")
new_react.name = "biomass with less metabolic requirement"
new_react.lower_bound = 0
new_react.upper_bound = 1000
new_react.reversibility = False
new_react_mets = dict()
for met in desired_biomass_components:
    new_react_met = model.metabolites.get_by_id(met)
    new_react_mets[new_react_met] = model.reactions.get_by_id(args.biomass_id).metabolites[model.metabolites.get_by_id(met)]
new_react.add_metabolites(new_react_mets)
reactions_toadd.add(new_react)
model.add_reactions(reactions_toadd)

# Remove old biomass reaction
model.remove_reactions(cobra.core.DictList([model.reactions.get_by_id(args.biomass_id)]))

# Set the new biomass reaction as objective
model.objective = "biomass_simplified"

# Calculate optimal growth with new biomass function
max_growth = model.slim_optimize()

# Obtain a minimal medium able to support at least 10% of the maximal biomass growth rate
required_medium = cobra.medium.minimal_medium(model, max_growth/10)

# Close all the medium's influxes
for reac_id in model.medium.keys():
    model.reactions.get_by_id(reac_id).lower_bound = 0.0

# Include the medium components specified by the user
for reac_id in desired_medium.keys():
    required_medium[reac_id] = desired_medium[reac_id]

# Set the new medium
for reac_id in required_medium.keys():
    model.reactions.get_by_id(reac_id).lower_bound = - required_medium[reac_id]

# Use fva to calculate the excretion and uptakes that the model can achieve while keeping 10% optimality
fva_sol = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0.1)
active_fva_ids = list(set(fva_sol[fva_sol["minimum"] < -1e-06].index)
                      .union(set(fva_sol[fva_sol["maximum"] > 1e-06].index)))
active_fva_exchanges = sorted([r_id for r_id in active_fva_ids if "EX_" in r_id])
possible_excretions = sorted([r_id for r_id in active_fva_exchanges if fva_sol.loc[r_id].maximum > 1e-06])
possible_uptakes = sorted([r_id for r_id in active_fva_exchanges if fva_sol.loc[r_id].minimum < -1e-06])

# Select all the sinks and demands that can carry a flux
boundary_ids = {r.id for r in model.boundary}
other_reacs_like_sinks_or_demands = sorted([r_id for r_id in boundary_ids.intersection(set(active_fva_ids)) 
                                              if r_id not in set(possible_excretions).union(set(possible_uptakes))])

# Remove from the model all the boundary reactions that are not feasible uptakes, excretions, sinks or demands
reacs_to_remove = list()
for r in model.boundary:
    if r.id not in set(possible_excretions).union(set(possible_uptakes)).union(other_reacs_like_sinks_or_demands):
        reacs_to_remove.append(r)        
model.remove_reactions(reacs_to_remove)

# Obtain a small list of reactions that support a 10% maximal biomass production(not proven to be minimal, but with 5 iterations, and using pFBA,the total number of reactions should be small). The excretion of all the compounds required by the user will be guaranteed. All medium reactions required by the user will be kept
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

# Set new id and name to the model
model.id = args.new_model_id
model.name = args.new_model_name

# Save the model
for reac in model.reactions:
    reac.lower_bound = float(reac.lower_bound)
    reac.upper_bound = float(reac.upper_bound)
cobra.io.write_sbml_model(model, args.output_file)
