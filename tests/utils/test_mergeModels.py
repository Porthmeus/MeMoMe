# SteFlor
# 15.05.24

import unittest

import cobra.io
import pandas as pd
from src.MeMoModel import *
from pathlib import Path
from src.mergeModels import *


class TestMergeModels(unittest.TestCase):

    def test_step1(self):
        tiny_model1 = "../dat/my_model_with_annotations.xml"
        tiny_model2 = "../dat/contiguous_pathway_model.xml"
        # read the input models
        model1 = MeMoModel.fromPath(Path(tiny_model1))
        model2 = MeMoModel.fromPath(Path(tiny_model2))
        matches = model1.match(model2)
        mo_me = ModelMerger(model1, model2, matches)
        cobra.io.write_sbml_model(mo_me.model1.cobra_model, "/home/stefano/Desktop/MeMoMe/tests/test_results/preproc_cobramodel1.xml")

        res = []
        self.assertTrue(all([res]))


if __name__ == '__main__':
    unittest.main()