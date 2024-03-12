# Porthmeus
# 24.02.23

import unittest
from pathlib import Path
import pandas as pd

from src.matchMets import matchMetsByInchi
from src.annotateInchiRoutines import inchiToMol, molToRDK, molToNormalizedInchi


class TestMatchInchit(unittest.TestCase):
    # The directory of this file
    this_directory = Path(__file__).parent
    dat = this_directory.joinpath("dat")

    def test_matchExpectedCases(self):
        # test the matchInchi algorithm and find expected matches
        # read the data of the test cases
        test_dat = pd.read_csv(self.dat.joinpath("InchiTestCases.csv"))
        res =[]
        for i in range(len(test_dat)):
                inchi1 = test_dat.loc[i,"Inchi1"]
                inchi2 = test_dat.loc[i,"Inchi2"]
                m1 = inchiToMol(inchi1)
                m2 = inchiToMol(inchi2)

                nminchi1 = molToNormalizedInchi(m1)
                nminchi2 = molToNormalizedInchi(m2)

                res.append(matchMetsByInchi(nminchi1, nminchi2, m1, m2, molToRDK(m1), molToRDK(m2))[0])
        expected_res = [bool(x) for x in test_dat.loc[:,"Result"]]
        self.assertTrue(all([x==y for x,y in zip(expected_res,res)]))

if __name__ == '__main__':
    unittest.main()
