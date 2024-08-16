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
        exchanges = model2.cobra_model.exchanges
        merger.convert_exchange_to_translation()
        self.assertTrue("t" in model2.cobra_model.compartments.keys())
        for ex in exchanges:
            self.assertTrue(ex.id.startswith("TR_"))
            self.assertEquals(len(list(ex.metabolites.keys())), 2)
            self.assertEquals({coef for coef in ex.metabolites.values()}, {-1.0, 1.0})
        merger.translate_reactions_and_metabolites_ids(score_thr=0.8)
        for ex in exchanges:
            print(ex.metabolites)



