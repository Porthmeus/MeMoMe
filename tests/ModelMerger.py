import unittest
from contextlib import redirect_stdout
from pathlib import Path
from copy import deepcopy

import cobra
import pandas as pd

from src.ModelMerger import ModelMerger
from src.MeMoModel import MeMoModel
from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix
from src.removeDuplicateMetabolites import removeDuplicateMetabolites

class Test_ModelMerger(unittest.TestCase):

    def setup(self):
        cmod1 = cobra.io.load_model("textbook")
        mod1 = MeMoModel().fromModel(cmod1)
    
        cmod2 = cobra.io.load_model("textbook")
        for met in cmod2.metabolites:
            met.id = "test2_" + met.id.replace("_e","_e0")
            if met.compartment == "e":
                met.compartment = "e0"
        for rxn in cmod2.exchanges:
            rxn.id =  rxn.id.replace("EX_","EX_test2_").replace("_e","_e0")
        cmod2._compartments.update({"e0":"external2"})

        cmod3 = cobra.io.load_model("textbook")
        for met in cmod3.metabolites:
            met.id = "test2_" + met.id.replace("_e","_e0")
            if met.compartment == "e":
                met.compartment = "e0"
        for rxn in cmod3.exchanges:
            rxn.id =  rxn.id.replace("EX_","EX_test3_").replace("_e","_e0")
        cmod2._compartments.update({"e0":"external3"})

        match_table = pd.DataFrame({"met_id1":sorted(list(set([handle_metabolites_prefix_suffix(x.id) for x in cmod1.metabolites]))),
                                    "met_id2":sorted(list(set([handle_metabolites_prefix_suffix(x.id) for x in cmod2.metabolites])))})
        return(match_table, cmod1, cmod2, cmod3)
    
    def test_translate(self):
        match_table, cmod1, cmod2, cmod3 = self.setup()
        # merge models
        mod1 = MeMoModel().fromModel(cmod1)
        mod2 = MeMoModel().fromModel(cmod2)
        mod1.annotated = True
        mod2.annotated = True
        modMerger = ModelMerger(mod1, mod2, match_table)
        modMerger.preprocess_models()
        modMerger.translate()
        modMerger.merge_models()

        self.assertEqual(modMerger.model1.compartments, {"c":"cytosol", "e":"extracellular", "t":"translation"})
        self.assertTrue(all([rxn.compartments == {"e"} for rxn in modMerger.model1.exchanges]))
        self.assertTrue(all([rxn.compartments == {"t","e"} for rxn in modMerger.model1.reactions if rxn.id.startswith("TR_")]))

        self.assertEqual(modMerger.model2.compartments, {"c":"cytosol", "e":"extracellular", "t":"translation"})
        self.assertTrue(all([rxn.compartments == {"e"} for rxn in modMerger.model2.exchanges]))
        self.assertTrue(all([rxn.compartments == {"t","e"} for rxn in modMerger.model2.reactions if rxn.id.startswith("TR_")]))
        
        # test whether the optimization result stays the same
        mergedMod = modMerger.merged_model.copy()
        for x in mergedMod.reactions:
            if x.lower_bound > 0:
                x.lower_bound = 0
        sol = mergedMod.optimize()
        for x in cmod1.reactions:
            if x.lower_bound > 0:
                x.lower_bound = 0
        for x in cmod2.reactions:
            if x.lower_bound > 0:
                x.lower_bound = 0
        self.assertTrue(round(sol.objective_value,5) == round(cmod1.optimize().objective_value,5))
        self.assertTrue(round(sol.objective_value,5) == round(cmod2.optimize().objective_value,5))
        self.assertEqual(mergedMod.medium, cmod1.medium)

        # test the prefix note in the merged model
        self.assertEqual(mergedMod.notes["MeMoMe_prefixes"],  {'M1': 'e_coli_core', 'M2': 'e_coli_core'})
        
        # test splitting of the models into single models again
        split_models = modMerger.split_merged_model()
        sol = split_models[0].optimize()
        ecore = cobra.io.load_model("textbook")
        sol_ecore = ecore.optimize()
        self.assertTrue(round(sol.objective_value,5) == round(sol_ecore.objective_value,5))
        sol = split_models[1].optimize()
        self.assertTrue(round(sol.objective_value,5) == round(sol_ecore.objective_value,5))
        self.assertTrue(all([True for x in split_models[0].metabolites if "e" in x.compartment and x.id[:-2] in list(match_table["met_id1"])]))
        self.assertTrue(all([True for x in split_models[1].metabolites if "e" in x.compartment and x.id[:-2] in list(match_table["met_id1"])]))

    def test_correct_prefix(self):

        match_table, cmod1, cmod2, cmod3 = self.setup()
        # merge models
        mod1 = MeMoModel().fromModel(cmod1)
        mod2 = MeMoModel().fromModel(cmod2)
        mod1.annotated = True
        mod2.annotated = True
        modMerger = ModelMerger(mod1, mod2, match_table)
        modMerger.preprocess_models()
        modMerger.translate()
        modMerger.merge_models()
        newMod1 = modMerger.merged_model
        # test the prefix translation
        self.assertTrue(any([x.id.startswith("M1") for x in newMod1.reactions]))
        newMod2 = modMerger.correct_prefix(newMod1, translation = {"M1":"M2","M2":"M3"})
        self.assertTrue(not any([x.id.startswith("M1") for x in newMod2.reactions]))
        self.assertTrue(not any([x.id.startswith("M1") for x in newMod2.metabolites]))
        self.assertTrue(not any([x.id.startswith("M2M1") for x in newMod2.reactions]))
        self.assertTrue(not any([x.id.startswith("M2M1") for x in newMod2.metabolites]))
        self.assertTrue(not any([x.id.startswith("M3M2") for x in newMod2.reactions]))
        self.assertTrue(not any([x.id.startswith("M3M2") for x in newMod2.metabolites]))
        
        # test merging of previous merged models
        mod1_exchangeMets = sorted(list(set([handle_metabolites_prefix_suffix(x.id) for x in newMod2.metabolites if x.compartment == "e"])))
        mod2_exchangeMets = ["test2_"+x for x in mod1_exchangeMets]
        match_table = pd.DataFrame({"met_id1":mod1_exchangeMets,

                                    "met_id2":mod2_exchangeMets})

        mod1 = MeMoModel().fromModel(newMod2)
        mod2 = MeMoModel().fromModel(cmod3)
        mod1.annotated = True
        mod2.annotated = True
        modMerger2 = ModelMerger(mod1, mod2, match_table)
        modMerger2.translate()
        modMerger2.merge_models()

        self.assertEqual(modMerger2.model1.compartments, {"c":"cytosol", "e":"extracellular", "t":"translation"})
        self.assertTrue(all([rxn.compartments == {"e"} for rxn in modMerger2.model1.exchanges]))
        self.assertTrue(all([rxn.compartments == {"t","e"} for rxn in modMerger2.model1.reactions if rxn.id.startswith("TR_")]))

        self.assertEqual(modMerger2.model2.compartments, {"c":"cytosol", "e":"extracellular", "t":"translation"})
        self.assertTrue(all([rxn.compartments == {"e"} for rxn in modMerger2.model2.exchanges]))
        self.assertTrue(all([rxn.compartments == {"t","e"} for rxn in modMerger2.model2.reactions if rxn.id.startswith("TR_")]))
        
        # test whether the optimization result stays the same
        mergedMod = modMerger2.merged_model.copy()
       # print(mergedMod.exchanges)
        self.assertEqual(mergedMod.medium, cmod1.medium)
        for x in mergedMod.reactions:
            if x.lower_bound > 0:
                x.lower_bound = 0
        sol = mergedMod.optimize()
        for x in cmod1.reactions:
            if x.lower_bound > 0:
                x.lower_bound = 0
        for x in cmod2.reactions:
            if x.lower_bound > 0:
                x.lower_bound = 0
        for x in cmod3.reactions:
            if x.lower_bound > 0:
                x.lower_bound = 0
        self.assertTrue(round(sol.objective_value,5) == round(cmod1.optimize().objective_value,5))
        self.assertTrue(round(sol.objective_value,5) == round(cmod2.optimize().objective_value,5))
        self.assertTrue(round(sol.objective_value,5) == round(cmod3.optimize().objective_value,5))
    
        # test the prefix note in the merged model
        self.assertEqual(mergedMod.notes["MeMoMe_prefixes"],  {'M1': 'e_coli_core', 'M2': 'e_coli_core', 'M3': 'e_coli_core'})


