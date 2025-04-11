import unittest
from src.MeMoModel import MeMoModel
from src.ModelMerger import ModelMerger
from pathlib import Path
import pandas as pd

class TestModelMerging(unittest.TestCase):
    this_directory = Path(__file__).parent
    dat = this_directory

    def test_translate_namespace(self):
        model1_path = self.dat.joinpath("dat/tiny_ecoli_keep_inchi.xml")
        model2_path = self.dat.joinpath("manually_merged_models/gapseq_recon3D/M2_bacterial_model.xml")
        model1 = MeMoModel.fromPath(model1_path)
        model2 = MeMoModel.fromPath(model2_path)
        #matches = model1.match(model2)
        # matches is for now hardcoded to simplify the test
        matches = pd.DataFrame({
            'met_id1': ['ac', 'ala_L'],
            'met_id2': ['cpd00029', 'cpd00035']
        })
        merger = ModelMerger(model1, model2, matches, fail_on_missing_metabolite=False)
        merger.translate_namespace()



