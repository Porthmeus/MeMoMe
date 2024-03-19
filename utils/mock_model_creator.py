import cobra
import random

# Function to generate a random InChI string (purely illustrative, not real InChI codes)
def generate_random_inchi():
    return f"InChI=1S/C{random.randint(1,20)}H{random.randint(1,40)}O{random.randint(1,10)}/p+1"

# Initialize the model
model = cobra.Model("example_model")

# Define metabolites
m1_e = cobra.Metabolite('M1_e', name='Metabolite1_ext', compartment='e')
m2_e = cobra.Metabolite('M2_e', name='Metabolite2_ext', compartment='e')
m3_e = cobra.Metabolite('M3_e', name='Metabolite3_ext', compartment='e')
m1_c = cobra.Metabolite('M1_c', name='Metabolite1', compartment='c')
m2_c = cobra.Metabolite('M2_c', name='Metabolite2', compartment='c')
m3_c = cobra.Metabolite('M3_c', name='Metabolite3', compartment='c')
m4_c = cobra.Metabolite('M4_c', name='Metabolite4', compartment='c')
m5_c = cobra.Metabolite('M5_c', name='Metabolite5', compartment='c')
m6_c = cobra.Metabolite('M6_c', name='Metabolite6', compartment='c')
m7_c = cobra.Metabolite('M7_c', name='Metabolite7', compartment='c')

# Define reactions operating on these metabolites

# Cytosolic reaction
reaction1 = cobra.Reaction('R1', lower_bound=-1000.0, upper_bound=1000.0)
reaction1.add_metabolites({m1_c: -1.0, m2_c: 2.0})

# Cytosolic reaction
reaction2 = cobra.Reaction('R2', lower_bound=-1000.0, upper_bound=1000.0)
reaction2.add_metabolites({m2_c: -2.0, m3_c: 3.0})

# Transport reaction
reaction3 = cobra.Reaction('R3', lower_bound=-1000.0, upper_bound=1000.0)
reaction3.add_metabolites({m1_e: -1.0, m1_c: 1.0}) 

# R4: Another cytosolic reaction
reaction4 = cobra.Reaction('R4', lower_bound=-1000.0, upper_bound=1000.0)
reaction4.add_metabolites({m2_c: -1.0, m4_c: 1.0})

# R5: A reaction producing m5_c from m4_c
reaction5 = cobra.Reaction('R5', lower_bound=-1000.0, upper_bound=1000.0)
reaction5.add_metabolites({m4_c: -2.0, m5_c: 2.0})

# R6: Metabolite conversion in cytosol
reaction6 = cobra.Reaction('R6', lower_bound=-1000.0, upper_bound=1000.0)
reaction6.add_metabolites({m5_c: -3.0, m6_c: 1.0})

# R7: Transport of m2 from external to cytosol
reaction7 = cobra.Reaction('R7', lower_bound=-1000.0, upper_bound=1000.0)
reaction7.add_metabolites({m2_e: -1.0, m2_c: 1.0})

# R8: External interaction, conversion of m3
reaction8 = cobra.Reaction('R8', lower_bound=-1000.0, upper_bound=1000.0)
reaction8.add_metabolites({m3_e: -2.0, m3_c: 2.0})

# R9: Synthesis of a new metabolite in cytosol from m6 and m7
reaction9 = cobra.Reaction('R9', lower_bound=-1000.0, upper_bound=1000.0)
reaction9.add_metabolites({m6_c: -1.0, m7_c: -1.0, m1_c: 1.0})

# R10: Export of m3 from cytosol to external
reaction10 = cobra.Reaction('R10', lower_bound=-1000.0, upper_bound=1000.0)
reaction10.add_metabolites({m3_c: -1.0, m3_e: 1.0})

# Add all new reactions to the model
model.add_reactions([reaction1, reaction2, reaction3, reaction4, reaction5,
                     reaction6, reaction7, reaction8, reaction9, reaction10])

# List of external metabolites for which to create exchange reactions
external_metabolites = [m1_e, m2_e, m3_e]

# Create and add exchange reactions for each external metabolite in a compact way
for metabolite in external_metabolites:
    exchange_reaction = cobra.Reaction(f'EX_{metabolite.id}', lower_bound=-1000.0, upper_bound=1000.0)
    exchange_reaction.add_metabolites({metabolite: -1.0})
    model.add_reactions([exchange_reaction])

# Adding annotations with random Seed IDs and InChI strings to each metabolite
for metabolite in model.metabolites:
    seed_id = f"cpd{random.randint(10000, 99999)}"
    inchi_string = generate_random_inchi()
    
    metabolite.annotation['seed.compound'] = seed_id
    metabolite.annotation['InChI'] = inchi_string

cobra.io.write_sbml_model(model, "../tests/dat/mock_model.xml")
