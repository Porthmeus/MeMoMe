# Porthmeus
# 24.02.23

import unittest
from pathlib import Path
from unittest import TestCase

import pandas as pd

from src.parseMetaboliteInfoFromSBML import *


class TestParseSBML(unittest.TestCase):
    # The directory of this file
    this_directory = Path(__file__).parent
    dat = this_directory.joinpath("dat")

    def test_parseExpectedValues(self):
        # test if one parses the metabolites from the ecoli core one gets the same metabolite table as expected
        # read the expected data frame
        ecoli_core_mets_df = self.dat.joinpath("ecoli_core_mets.csv")
        expected_df = pd.read_csv(ecoli_core_mets_df, index_col = 0)
        
        # parse the model
        ecoli_core_sbml = self.dat.joinpath("e_coli_core.xml")
        test_df = parseMetaboliteInfoFromSBML(ecoli_core_sbml)
        self.assertTrue(all(test_df == expected_df))
    
    def test_SBMLfileNonExistant(self):
        # test whether the function still handles non existant files, as libSBML will not complain about that by default
        non_existent_path = Path("nonexsistent")
        self.assertRaises(FileExistsError, parseMetaboliteInfoFromSBML,non_existent_path)


    def test_extract_cpd(self):
        teststr = "M_cpd00001_c0"
        self.assertEqual(extract_cpd(teststr), "cpd00001")
        teststr = "M_cPd00001_c0"
        self.assertEqual(extract_cpd(teststr), None)


if __name__ == '__main__':
    unittest.main()
