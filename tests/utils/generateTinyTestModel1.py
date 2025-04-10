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

# Internal pyruvate
pyr_c = Metabolite('pyr_c', formula='C3H3O3', name='Pyruvate', compartment='c')
pyr_c.annotation['bigg.metabolite'] = 'pyr'

# External pyruvate
pyr_e = Metabolite('pyr_e', formula='C3H3O3', name='Pyruvate', compartment='e')
pyr_e.annotation['bigg.metabolite'] = 'pyr'

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

# Transport of pyruvate from intracellular to extracellular
rxn3 = Reaction('T_pyr')
rxn3.name = 'Pyruvate Transport'
rxn3.lower_bound = 0
rxn3.upper_bound = 1000
rxn3.add_metabolites({
    pyr_c: -1.0,  # Internal pyruvate consumption
    pyr_e: 1.0    # External pyruvate production
})

# Define exchange reactions for external metabolites
# ATP exchange
atp_exchange = Reaction('EX_atp_e')
atp_exchange.name = 'ATP exchange'
atp_exchange.lower_bound = -10  # Arbitrary lower bound for uptake
atp_exchange.upper_bound = 1000
atp_exchange.add_metabolites({
    atp_e: -1.0  # ATP uptake from the extracellular space
})

# ADP exchange
adp_exchange = Reaction('EX_adp_e')
adp_exchange.name = 'ADP exchange'
adp_exchange.lower_bound = -10  # Arbitrary lower bound for uptake
adp_exchange.upper_bound = 1000
adp_exchange.add_metabolites({
    adp_e: -1.0  # ADP uptake from the extracellular space
})

# Pyruvate exchange
pyr_exchange = Reaction('EX_pyr_e')
pyr_exchange.name = 'Pyruvate exchange'
pyr_exchange.lower_bound = -10  # Arbitrary lower bound for uptake
pyr_exchange.upper_bound = 1000
pyr_exchange.add_metabolites({
    pyr_e: -1.0  # Pyruvate uptake from the extracellular space
})

# Define reactions within the model
# Reaction 1: Hexokinase (glucose to glucose-6-phosphate)
glucose = Metabolite('glc__D_c', formula='C6H12O6', name='D-Glucose', compartment='c')
glucose.annotation['bigg.metabolite'] = 'glc__D'
g6p_c = Metabolite('g6p_c', formula='C6H13O9P', name='Glucose-6-phosphate', compartment='c')
g6p_c.annotation['bigg.metabolite'] = 'g6p'

hexokinase = Reaction('hexokinase')
hexokinase.name = 'Hexokinase'
hexokinase.lower_bound = 0
hexokinase.upper_bound = 1000
hexokinase.add_metabolites({
    glucose: -1.0,  # Glucose consumption
    atp_c: -1.0,    # ATP consumption
    g6p_c: 1.0,     # Glucose-6-phosphate production
    adp_c: 1.0      # ADP production
})

# Reaction 2: Pyruvate production from glucose-6-phosphate
glycolysis = Reaction('glycolysis')
glycolysis.name = 'Glycolysis'
glycolysis.lower_bound = 0
glycolysis.upper_bound = 1000
glycolysis.add_metabolites({
    g6p_c: -1.0,  # Glucose-6-phosphate consumption
    pyr_c: 1.0    # Pyruvate production
})

# Add all metabolites and reactions to the model
model.add_metabolites([atp_c, atp_e, adp_c, adp_e, pyr_c, pyr_e, glucose, g6p_c])
model.add_reactions([rxn1, rxn2, rxn3, hexokinase, glycolysis, atp_exchange, adp_exchange, pyr_exchange])

# Print model summary
print(model.summary())

# Save the model as XML
io.write_sbml_model(model, "tests/dat/my_model_with_annotations.xml")
