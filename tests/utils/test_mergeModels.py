# SteFlor
# 15.05.24

import unittest

import cobra.io
import pandas as pd
from pathlib import Path
import src.mergeModels

class TestmergeModels(unittest.TestCase):

    def test_prepend_ids(self):
        #read the input model
        model = cobra.io.read_sbml_model("tests/dat/my_model_with_annotations.xml")

        res = []
        self.assertTrue(all([res]))




if __name__ == '__main__':
    unittest.main()