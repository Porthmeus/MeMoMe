#Hadjar
#Created tests
#10.05.2023

# Stefano
# Allowed automatic testing
# 11.05.2023

import unittest
from pathlib import Path

from src.removeCompartment import removeCompAndOrganismFromMetIds



class TestRemoveCompartment(unittest.TestCase):

    def test_matchExpectedCases(self):
        # test the matchInchi algorithm and find expected matches
        # define the test cases and respective results
        metabolite_ids = ['M_cpd00001__91__c__93__', 'M_cpd00002__lparen__e__rparen__', 'M_cpd00003__91__c__93__','M_cpd00001_c', 'M_cpd00002_e', 'M_cpd00003_p', 'M_cpd00004_c','glucose[c]', 'oxygen[e]', 'atp[c]']
        expected_new_metabolite_ids = ['cpd00001', 'cpd00002', 'cpd00003','cpd00001', 'cpd00002', 'cpd00003', 'cpd00004','glucose', 'oxygen', 'atp']

        res = removeCompAndOrganismFromMetIds(metabolite_ids)
        self.assertTrue(all([x==y for x,y in zip(expected_new_metabolite_ids, res)]))

if __name__ == '__main__':
    unittest.main()
