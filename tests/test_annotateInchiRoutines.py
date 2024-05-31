from unittest import TestCase

from src.annotateInchiRoutines import findOptimalInchi


class TestInchiRoutines(TestCase):
    def test_find_optimal_inchi(self):
        # test invalid inchis
        s = ['InChI=1S/C5H13NO4S/c1-6(2,3)4-5-10-11(7,8)9/h4-5H2,1-3H3',
             'InChI=1S/C5H13NO4S/c1-6(2,3)4-5-10-11(7,8)9/h4-5H2,1-3H3/p+1']
        self.assertTrue(findOptimalInchi(s) is None)
        
        # test whether the longest inchi is correctly returned and ignore charge rule, if charge == none
        s = ["InChI=1S/C5H11O8P/c6-3-2(1-12-14(9,10)11)13-5(8)4(3)7/h2-8H,1H2,(H2,9,10,11)/t2-,3-,4-,5+/m1/s1",
             "InChI=1S/C5H11O8P/c6-3-2(1-12-14(9,10)11)13-5(8)4(3)7/h2-8H,1H2,(H2,9,10,11)/p-2/t2-,3-,4-,5+/m1/s1"]
        self.assertTrue(findOptimalInchi(s, charge = None) == s[1])

        # test whether charge rule works
        self.assertTrue(findOptimalInchi(s, charge = 0) == s[0])

        # test whether charge is ignored if none of the inchis has the correct charge
        self.assertTrue(findOptimalInchi(s, charge = 8) == s[1])

