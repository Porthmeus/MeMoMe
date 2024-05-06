import cobra
from cobra import Model, Metabolite, Reaction, io

# Create a new COBRApy model
model = Model('my_model_with_annotations')

# Define internal and external metabolites
# Internal ATP
atp_c = Metabolite('atp_c', formula='C10H16N5O13P3', name='ATP', compartment='c')
atp_c.annotation['bigg.metabolite'] = 'atp'

# External ATP
atp_e = Metabolite('atp_e', formula='C10H16N5O13P3', name='ATP', compartment='e')
atp_e.annotation['bigg.metabolite'] = 'atp'

# Internal ADP
adp_c = Metabolite('adp_c', formula='C10H15N5O10P2', name='ADP', compartment='c')
adp_c.annotation['bigg.metabolite'] = 'adp'

# External ADP
adp_e = Metabolite('adp_e', formula='C10H15N5O10P2', name='ADP', compartment='e')
adp_e.annotation['bigg.metabolite'] = 'adp'

# Internal glucose-6-phosphate
g6p_c = Metabolite('g6p_c', formula='C6H13O9P', name='Glucose-6-phosphate', compartment='c')
g6p_c.annotation['bigg.metabolite'] = 'g6p'

# External glucose-6-phosphate
g6p_e = Metabolite('g6p_e', formula='C6H13O9P', name='Glucose-6-phosphate', compartment='e')
g6p_e.annotation['bigg.metabolite'] = 'g6p'

# Define transport reactions between compartments
# Transport of ATP from extracellular to intracellular
rxn1 = Reaction('T_atp')
rxn1.name = 'ATP Transport'
rxn1.lower_bound = 0
rxn1.upper_bound = 1000
rxn1.add_metabolites({
    atp_e: -1.0,  # External ATP consumption
    atp_c: 1.0    # Internal ATP production
})

# Transport of ADP from intracellular to extracellular
rxn2 = Reaction('T_adp')
rxn2.name = 'ADP Transport'
rxn2.lower_bound = 0
rxn2.upper_bound = 1000
rxn2.add_metabolites({
    adp_c: -1.0,  # Internal ADP consumption
    adp_e: 1.0    # External ADP production
})

# Transport of glucose-6-phosphate from extracellular to intracellular
rxn3 = Reaction('T_g6p')
rxn3.name = 'Glucose-6-phosphate Transport'
rxn3.lower_bound = 0
rxn3.upper_bound = 1000
rxn3.add_metabolites({
    g6p_e: -1.0,  # External glucose-6-phosphate consumption
    g6p_c: 1.0    # Internal glucose-6-phosphate production
})

# Add the metabolites and reactions to the model
model.add_metabolites([atp_c, atp_e, adp_c, adp_e, g6p_c, g6p_e])
model.add_reactions([rxn1, rxn2, rxn3])

# Print model summary
print(model.summary())

# Save the model to a xml file
io.write_sbml_model(model,"tests/dat/my_model_with_annotations.xml")