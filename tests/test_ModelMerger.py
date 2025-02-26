import unittest
from src.MeMoModel import MeMoModel
from src.ModelMerger import ModelMerger
from pathlib import Path
import pandas as pd
import cobra
import re


def create_minimal_model(model_id: str, id_met1: str, id_met2: str):
    """
    Create a minimal COBRA model with two metabolites and basic transport/exchange reactions.

    Parameters:
        model_id (str): Identifier for the COBRA model.
        id_met1 (str): Identifier for the first metabolite.
        id_met2 (str): Identifier for the second metabolite.

    Returns:
        cobra.Model: A COBRA model with minimal reactions.
    """
    model = cobra.Model(model_id)

    # Define metabolites
    # Cytosolic metabolites
    met1_c = cobra.Metabolite(id_met1 + "_c", compartment="c")
    met2_c = cobra.Metabolite(id_met2 + "_c", compartment="c")
    # External metabolite
    met1_e = cobra.Metabolite(id_met1 + "_e", compartment="e")

    # Define reactions
    # Conversion reaction: Converts metabolite 1 to metabolite 2 in the cytosol
    r1 = cobra.Reaction("R1")
    r1.name = "Conversion of " + id_met1 + " to " + id_met2
    r1.add_metabolites({met1_c: -1, met2_c: 1})
    r1.bounds = (-1000, 1000)  # Allow bidirectional reaction

    # Transport reaction: Moves metabolite 1 from cytosol to external compartment
    r2 = cobra.Reaction("R2")
    r2.name = "Transport of " + id_met1 + " from cytosol to external"
    r2.add_metabolites({met1_c: -1, met1_e: 1})
    r2.bounds = (-1000, 1000)  # Allow bidirectional transport

    # Exchange reaction: Represents uptake/secretion of metabolite 1 in the external medium
    ex = cobra.Reaction("EX_" + id_met1 + "_e")
    ex.name = id_met1 + " exchange"
    ex.add_metabolites({met1_e: -1})  # Consumption of the external metabolite

    # Add reactions to the model
    model.add_reactions([r1, r2, ex])

    return model


class TestModelMerging(unittest.TestCase):
    this_directory = Path(__file__).parent
    dat = this_directory

    # create model1 using cobrapy
    model1 = create_minimal_model("model1","a","b")
    model2 = create_minimal_model("model2","c","d")

    # save models as temporary xml files
    cobra.io.write_sbml_model(model1, "utils/tmp/model1.xml")
    cobra.io.write_sbml_model(model2, "utils/tmp/model2.xml")


    def test_automatic_compartment_creation(self):
        """
        Makes sure that the assignment of a non-existing compartment to a metabolite results in the creation of the compartment
        """
        model_path = self.dat.joinpath("dat/tiny_ecoli_keep_inchi.xml")
        model = MeMoModel.fromPath(model_path)
        cobra_model = model.cobra_model
        # ensure that the model that I am using to test this cobrapy functionality doesn't already have a "t" compartment
        assert "t" not in cobra_model.compartments.keys()
        # assign a metabolite to the translation compartment. This compartment still does not exist, but it will be
        # created automatically by assigning a metabolite to it
        met = cobra_model.metabolites[0]
        met.compartment = "t"
        # ensure that the translation compartment has been added
        assert "t" in cobra_model.compartments.keys()


    def test_translate_namespace(self):
        model1_path = self.dat.joinpath("dat/tiny_ecoli_keep_inchi.xml")
        model2_path = self.dat.joinpath("manually_merged_models/gapseq_recon3D/M2_bacterial_model.xml")
        model1 = MeMoModel.fromPath(model1_path)
        model2 = MeMoModel.fromPath(model2_path)
        matches = model1.match(model2)
        merger = ModelMerger(model2, matches)
        pre_translation_exchanges = model2.cobra_model.exchanges
        merger.convert_exchange_rxns_to_translation_rxns()
        self.assertTrue("t" in model2.cobra_model.compartments.keys())
        for ex in pre_translation_exchanges:
            self.assertTrue(ex.id.startswith("TR_"))
            self.assertEqual(len(list(ex.metabolites.keys())), 2)
            self.assertEqual({coef for coef in ex.metabolites.values()}, {-1.0, 1.0})
        merger.translate_rxn_and_met_ids(score_thr=0.5)
        # TODO: manually create a complete list of "correct matches" between the 2 tested models, and test if they are assigned correctly
        merger.create_exchanges()
        for pre_ex in pre_translation_exchanges:
            met_id = pre_ex.id.replace("TR_","")
            met = model2.cobra_model.metabolites.get_by_id(met_id)
            new_ex = model2.cobra_model.reactions.get_by_id("EX_" + met_id)
            self.assertEqual(type(new_ex.metabolites), dict)
            self.assertEqual(len(new_ex.metabolites), 1)
            self.assertEqual([coef for coef in new_ex.metabolites.values()], [-1.0])
            self.assertEqual([key for key in new_ex.metabolites.keys()], [met])
        bounds = dict()
        for rxn in pre_translation_exchanges:
            bounds[rxn] = (rxn.lower_bound, rxn.upper_bound)
        merger.set_rxn_bounds()
        for rxn in pre_translation_exchanges:
            self.assertEqual(rxn.lower_bound, -1000)
            self.assertEqual(rxn.upper_bound, 1000)
            new_ex = model2.cobra_model.reactions.get_by_id(rxn.id.replace("TR_", "EX_"))
            self.assertEqual(bounds[rxn], (new_ex.lower_bound, new_ex.upper_bound))


    def test_case_1(self):

        # load each of them into a MeMoModel with MeMoModel.fromPath(model_path)
        model_path1 = self.this_directory.joinpath("tmp/model1.xml")
        model_path2 = self.this_directory.joinpath("tmp/model2.xml")
        meMoModel1 = MeMoModel.fromPath(model_path1)
        meMoModel2 = MeMoModel.fromPath(model_path2)

        #define an empty matching table
        matches = pd.DataFrame()
        # function call
        merge(meMoModel1, meMoModel2, 1, 2, matches)

        # check if the set of reaction ids of the merged model correspond to what expected in the test
        expected_rxn_ids = set()

        # check if the set of metabolite ids of the merged model correspond to what expected in the test
        expected_met_ids = set()



