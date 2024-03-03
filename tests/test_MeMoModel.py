# Porthmeus
# 02.08.23

import unittest
from pathlib import Path
import cobra as cb
import pandas as pd
import warnings

#from src.MeMoModel import MeMoModel
from src.MeMoMetabolite import MeMoMetabolite
from src.MeMoModel import MeMoModel
from src.annotateModelSEED import annotateModelSeed, annotateModelSeed_id
from src.annotateChEBI import annotateChEBI
from src.annotateBiGG import annotateBiGG, annotateBiGG_id

class Test_annotateBulkRoutines(unittest.TestCase):
    # The directory of this file
    #this_directory = Path("tests")
    this_directory = Path(__file__).parent
    dat = this_directory.joinpath("dat")

    def test_MeMoModel(self):
        # load e.coli core as a reference in two of the three supported methods
        mod_path = self.dat.joinpath("e_coli_core.xml")
        mod = MeMoModel.fromPath(mod_path)
        self.assertIsInstance(mod, MeMoModel)
        cb_mod = cb.io.read_sbml_model(str(mod_path))
        mod = MeMoModel.fromModel(cb_mod)
        self.assertIsInstance(mod, MeMoModel)

    def test_MeMoModelAnnotation(self):
        # load the e.coli core model and bulk annotate the metabolites. Check if any annoation tkes place (Chebi should cover all metabolites)
        mod_path = self.dat.joinpath("e_coli_core.xml")
        mod = MeMoModel.fromPath(mod_path)
        pre_inchis = [x._inchi_string for x in mod.metabolites]
        mod.annotate()
        post_inchis = [x._inchi_string for x in mod.metabolites]
        self.assertTrue(any([x != y for x,y in zip(pre_inchis, post_inchis)]))
        # check if it also works if we remove the annotations
        mod = MeMoModel.fromPath(mod_path)
        # ignore the warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            [x.set_annotations({}) for x in mod.metabolites]
        mod.annotate()

    def test_annotateChEBI(self):
        # test the matchInchi algorithm and find expected matches
        # read the data of the test cases - I picked three random cases from the ChEBI file
        test_dat = ["2873", "5850", "176363"]
        inchis = ["InChI=1S/C30H48O5/c1-17-9-12-30(25(34)35)14-13-28(5)19(23(30)18(17)2)7-8-22-26(3)15-20(32)24(33)27(4,16-31)21(26)10-11-29(22,28)6/h7,17-18,20-24,31-33H,8-16H2,1-6H3,(H,34,35)/t17-,18+,20-,21-,22-,23+,24+,26+,27+,28-,29-,30+/m1/s1", "InChI=1S/C18H16N2O2/c1-22-17-9-5-4-8-16(17)20-18(21)19-15-11-10-13-6-2-3-7-14(13)12-15/h2-12H,1H3,(H2,19,20,21)","InChI=1S/C26H28O14/c27-6-12-17(31)21(35)23(37)25(39-12)15-19(33)14-10(30)5-11(8-1-3-9(29)4-2-8)38-24(14)16(20(15)34)26-22(36)18(32)13(7-28)40-26/h1-5,12-13,17-18,21-23,25-29,31-37H,6-7H2/t12-,13+,17-,18+,21+,22-,23-,25?,26?/m1/s1"]
        metabolites = []
        for chebi in test_dat:
            met = MeMoMetabolite()
            dic = {"chebi":[chebi]}
#            with warnings.catch_warnings(action = "ignore"):
                # TODO figure out why this is throwing an error!
            met.set_annotations(dic)
            metabolites.append(met)
        
        annotateChEBI(metabolites)
        self.assertTrue(all([y==z for y,z in zip([x._inchi_string for x in metabolites], inchis)]))    

    def test_annotateBiGG(self):
        # create a small test for the annotateBigg functions
        # create a mock list of metabolites
        m1 = MeMoMetabolite(_id = "glc__D")
        m2 = MeMoMetabolite(_id = "mock_id",annotations = {"bigg.metabolite":["glc__D"]})
        mets = [m1,m2]
        x,y,z = annotateBiGG(mets)
        self.assertTrue(x == 0)
        self.assertTrue(y == 1)
        self.assertTrue(z == 1)
        x,y,z = annotateBiGG_id(mets)
        self.assertTrue(x == 0)
        self.assertTrue(y == 1)
        self.assertTrue(z == 1)
        # redo to test for correct counting
        x,y,z = annotateBiGG(mets)
        self.assertTrue(x == 0)
        self.assertTrue(y == 0)
        self.assertTrue(z == 0)
        x,y,z = annotateBiGG_id(mets)
        self.assertTrue(x == 0)
        self.assertTrue(y == 0)
        self.assertTrue(z == 0)

    def test_annotateModelSEED(self):
        # create a small test for the annotateBigg functions
        # create a mock list of metabolites
        m1 = MeMoMetabolite(_id = "cpd00027")
        m2 = MeMoMetabolite(_id = "mock_id",annotations = {"seed.compound":["cpd00027"]})
        mets = [m1,m2]
        x,y,z = annotateModelSeed(mets)
        self.assertTrue(x == 1)
        self.assertTrue(y == 1)
        self.assertTrue(z == 1)
        x,y,z = annotateModelSeed_id(mets)
        self.assertTrue(x == 1)
        self.assertTrue(y == 1)
        self.assertTrue(z == 1)
        # redo the test and check that nothing is added
        x,y,z = annotateModelSeed_id(mets)
        self.assertTrue(x == 0)
        self.assertTrue(y == 0)
        self.assertTrue(z == 0)
        x,y,z = annotateModelSeed(mets)
        self.assertTrue(x == 0)
        self.assertTrue(y == 0)
        self.assertTrue(z == 0)

    def test_MeMoModelCompare(self):
        # test the comparison for metabolite matching
        mod_path = self.dat.joinpath("e_coli_core.xml")
        mod = MeMoModel.fromPath(mod_path)
        # self comparison
        res = mod.match(mod)
        self.assertIsInstance(res, pd.DataFrame)
        self.assertTrue(all([x in res.columns for x in ["met_id1","met_id2"]]))

if __name__ == '__main__':
    unittest.main()
