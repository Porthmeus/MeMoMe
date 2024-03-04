from unittest import TestCase

from src.annotateInchiRoutines import findOptimalInchi


class TestInchiRoutines(TestCase):
    def test_find_optimal_inchi(self):
        s = ['InChI=1S/C5H13NO4S/c1-6(2,3)4-5-10-11(7,8)9/h4-5H2,1-3H3',
             'InChI=1S/C5H13NO4S/c1-6(2,3)4-5-10-11(7,8)9/h4-5H2,1-3H3/p+1']
        self.assertTrue(findOptimalInchi(s) is None)
