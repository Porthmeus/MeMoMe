import unittest
from src.MeMoModel import MeMoModel
from src.ModelMerger import ModelMerger
from pathlib import Path
import re

class TestModelMerging(unittest.TestCase):
    this_directory = Path(__file__).parent
    dat = this_directory

    def test_translate_namespace(self):
        model1_path = self.dat.joinpath("dat/tiny_ecoli_keep_inchi.xml")
        model2_path = self.dat.joinpath("manually_merged_models/gapseq_recon3D/M2_bacterial_model.xml")
        model1 = MeMoModel.fromPath(model1_path)
        model2 = MeMoModel.fromPath(model2_path)
        matches = model1.match(model2)
        merger = ModelMerger(model2, matches)
        pre_translation_exchanges = model2.cobra_model.exchanges
        merger.replace_exchange_rxns_with_translation_rxns()
        self.assertTrue("t" in model2.cobra_model.compartments.keys())
        for ex in pre_translation_exchanges:
            self.assertTrue(ex.id.startswith("TR_"))
            self.assertEqual(len(list(ex.metabolites.keys())), 2)
            self.assertEqual({coef for coef in ex.metabolites.values()}, {-1.0, 1.0})
        merger.translate_ids()
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



