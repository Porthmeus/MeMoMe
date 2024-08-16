# Porthmeus
# 24.02.23
# coding: utf-8

import unittest
from pathlib import Path
from src.MeMoModel import MeMoModel
from src.MeMoMetabolite  import MeMoMetabolite



class Test_MeMoMetabolite(unittest.TestCase):

    def test_emptyInchiStringHandling(self):
        b = MeMoModel()
        metabolite: MeMoMetabolite = MeMoMetabolite()
        metabolite.set_id("A")
        metabolite.set_inchi_string("")
        metabolite.add_inchi_string("")
        metabolite.set_names(["Peter"])
        self.assertTrue(metabolite._inchi_string == None)
        

    def test_AnnotateMatchingOnEmptyInchis(self):
        b = MeMoModel()

        metabolite: MeMoMetabolite = MeMoMetabolite()
        metabolite.set_id("A")
        metabolite.set_inchi_string("")
        metabolite.add_inchi_string("")
        metabolite.set_names(["Peter"])


        metaboliteB: MeMoMetabolite = MeMoMetabolite()
        metaboliteB.set_id("B")
        metaboliteB.set_inchi_string("")
        metaboliteB.add_inchi_string("")
        metaboliteB.set_names(["Peter"])

        model1 = MeMoModel([metabolite])
        model2 = MeMoModel([metaboliteB])

        model1.annotate()
        model2.annotate()

        print(model1.match(model2))

if __name__ == "__main__":
    unittest.main()
