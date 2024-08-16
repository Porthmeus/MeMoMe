import unittest
from src.MeMoModel import MeMoModel
from src.ModelMerger import ModelMerger
from pathlib import Path

class TestModelMerging(unittest.TestCase):
    this_directory = Path(__file__).parent
    dat = this_directory

    def test_translate_namespace(self):
        model1_path = self.dat.joinpath("dat/tiny_ecoli_keep_inchi.xml")
        model2_path = self.dat.joinpath("manually_merged_models/gapseq_recon3D/M2_bacterial_model.xml")
        model1 = MeMoModel.fromPath(model1_path)
        model2 = MeMoModel.fromPath(model2_path)
        matches = model1.match(model2)
        print(matches)
        merger = ModelMerger(model2, matches)
        merger.convert_exchange_to_translation()
        print(model2.cobra_model.compartments)
        self.assertTrue("t" in model2.cobra_model.compartments.keys())
        for reac in model2.cobra_model.exchanges:
            self.assertFalse(reac.id.startswith("EX_"))
        merger.translate_reactions_and_metabolites_ids()
        for reac in model2.cobra_model.exchanges:
            print(reac.metabolites)


