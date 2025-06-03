# Porthmeus
# 22.05.25 

import unittest
from src.dbhandling.reformatAux import *
from pathlib import Path
import pandas as pd
import os
import shutil

class Test_DBreformatting(unittest.TestCase):
    # The directory of this file
    #this_directory = Path("tests")
    this_directory = Path(__file__).parent
    dat = this_directory.joinpath("dat")

    def setUp(self):
        self.test_dir = "TestCaseDir"
        os.makedirs(self.test_dir, exist_ok = True)

    def tearDown(self):
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def test_getData(self):
        # test that unallowed dbs lead to an error
        with self.assertRaises(ValueError):
            getData("foobar")
        
        # test if all implemented databases return pd.dataframes
        self.assertEqual(type(getData("VMH")), type(pd.DataFrame()))
        self.assertEqual(type(getData("BiGG")), type(pd.DataFrame()))
        self.assertEqual(type(getData("ModelSeed")), type(pd.DataFrame()))

    def test_writeData(self):
        # test if we can create an empty file
        writeData(pd.DataFrame(), db = "TestCase")
        self.assertTrue(os.path.exists(os.path.join(self.test_dir,"TestCaseFile.txt")))

