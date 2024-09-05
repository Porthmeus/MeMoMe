# Porthmeus
# 02.08.23

import unittest
import cobra as cb
import pandas as pd
import warnings
import math
import sys
from pathlib import Path
from src.MeMoMetabolite import MeMoMetabolite
from src.MeMoModel import MeMoModel
from src.annotateModelSEED import annotateModelSEED, annotateModelSEED_id
from src.annotateChEBI import annotateChEBI
from src.annotateBiGG import annotateBiGG, annotateBiGG_id
from src.annotateAux import AnnotationResult
from src.removeDuplicateMetabolites import detectDuplicates, removeDuplicateMetabolites

print(sys.version)

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
        self.assertIsInstance(mod.cobra_model, cb.Model)
        cb_mod = cb.io.read_sbml_model(str(mod_path))
        mod = MeMoModel.fromModel(cb_mod)
        self.assertIsInstance(mod, MeMoModel)
        self.assertIsInstance(mod.cobra_model, cb.Model)

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
        
        anno_res = annotateChEBI(metabolites)
        self.assertTrue(anno_res == AnnotationResult(3,0,0))
        self.assertTrue(all([y==z for y,z in zip([x._inchi_string for x in metabolites], inchis)]))    

    def test_annotateBiGG(self):
        # create a small test for the annotateBigg functions
        # create a mock list of metabolites
        m1 = MeMoMetabolite(_id = "glc__D")
        m2 = MeMoMetabolite(_id = "mock_id",annotations = {"bigg.metabolite":["glc__D"]})
        mets = [m1,m2]
        anno_res = annotateBiGG(mets)
        self.assertTrue(anno_res == AnnotationResult(0,1,1))
        anno_res = annotateBiGG_id(mets)
        self.assertTrue(anno_res == AnnotationResult(0,1,1))
        # redo to test for correct counting
        anno_res = annotateBiGG(mets)
        self.assertTrue(anno_res == AnnotationResult(0,0,0))
        anno_res = annotateBiGG_id(mets)
        self.assertTrue(anno_res == AnnotationResult(0,0,0))

    def test_annotateModelSEED(self):
        # create a small test for the annotateBigg functions
        # create a mock list of metabolites
        m1 = MeMoMetabolite(_id = "cpd00027")
        m2 = MeMoMetabolite(_id = "mock_id",annotations = {"seed.compound":["cpd00027"]})
        mets = [m1,m2]
        anno_res = annotateModelSEED(mets)
        self.assertTrue(anno_res == AnnotationResult(1,1,1))
        anno_res = annotateModelSEED_id(mets)
        self.assertTrue(anno_res == AnnotationResult(1,1,1))
        # redo the test and check that nothing is added
        anno_res = annotateModelSEED_id(mets)
        self.assertTrue(anno_res == AnnotationResult(0,0,0))
        anno_res = annotateModelSEED(mets)
        self.assertTrue(anno_res == AnnotationResult(0,0,0))

    def test_MeMoModelCompare(self):
        # test the comparison for metabolite matching
        mod_path = self.dat.joinpath("e_coli_core.xml")
        mod = MeMoModel.fromPath(mod_path)
        mod.annotate()
        # self comparison
        res = mod.match(mod,keepAllMatches = False)
        self.assertIsInstance(res, pd.DataFrame)
        self.assertTrue(all([x in res.columns for x in ["met_id1","met_id2"]]))
        self.assertTrue(all([x==y  for x,y in zip(res.met_id1,res.met_id2)]))


class Test_MiscStuff(unittest.TestCase):

    def test_MeMoModelOutPutNames(self):
      metaboliteA: MeMoMetabolite = MeMoMetabolite()
      metaboliteB: MeMoMetabolite = MeMoMetabolite()
      metaboliteA.set_names(["Glucose"])
      metaboliteB.set_names(["Glukose"])

      model = MeMoModel([metaboliteA])
      model2 = MeMoModel([metaboliteB])
      res = model.match(model2, output_names = False)
      self.assertFalse("name_id1" in res.columns)
      self.assertFalse("name_id2" in res.columns)
      
      res = model.match(model2, output_names = True)
      self.assertTrue("name_id1" in res.columns)
      self.assertTrue("name_id2" in res.columns)

      self.assertEqual(res["name_id1"][0], "Glucose")
      self.assertEqual(res["name_id2"][0], "Glukose")


    def test_MeMoModelOutputDBs(self):
      metaboliteA: MeMoMetabolite = MeMoMetabolite()
      metaboliteB: MeMoMetabolite = MeMoMetabolite()
      # TODO: How does a sensible annoation example look like 
      metaboliteA.set_annotations({"DatabaseA" : ["stuff", "stuff3"]})
      metaboliteB.set_annotations({"DatabaseA" : ["stuff"]})

      model = MeMoModel([metaboliteA])
      model2 = MeMoModel([metaboliteB])
      res = model.match(model2, output_names = False, output_dbs = False)
      self.assertFalse("commonDBs" in res.columns)
      self.assertFalse("commonIds" in res.columns)
      self.assertFalse("distinctIds" in res.columns)

      res = model.match(model2, output_names = False, output_dbs = True)
      self.assertTrue("commonDBs" in res.columns)
      self.assertTrue("commonIds" in res.columns)
      self.assertTrue("allIds" in res.columns)
      
      self.assertEqual(res["commonDBs"][0], "DatabaseA")
      self.assertEqual(res["commonIds"][0], "['DatabaseA.stuff']")
      self.assertEqual(res["allIds"][0], "['DatabaseA.stuff', 'DatabaseA.stuff3']")

    def test_1toManyMatchingOnName(self):
      metaboliteA: MeMoMetabolite = MeMoMetabolite()
      metaboliteB: MeMoMetabolite = MeMoMetabolite()
      metaboliteA.set_names(["Glucose"])
      metaboliteB.set_names(["Glukose"])
      metaboliteA.set_id("A")
      metaboliteB.set_id("B")

      model = MeMoModel([metaboliteA, metaboliteB])
      model2 = MeMoModel([metaboliteA, metaboliteB])
      res = model.match(model2, keepAllMatches = True)

      self.assertEqual(res.shape[0], 4)
      self.assertEqual(res["Name_score"][0], 1.0000)
      val = res["Name_score"][1]
      self.assertTrue(math.isclose(val, 0.857143, rel_tol=1e-2))

      res = model.match(model2, keepAllMatches = False)
      self.assertEqual(res.shape[0], 2)
      self.assertEqual(res["Name_score"][0], 1.0000)
      self.assertEqual(res["Name_score"][1], 1.0000)
      
      res = model.match(model2, keepAllMatches = True, threshold_name = 0.9)
      self.assertEqual(res.shape[0], 2)
      self.assertEqual(res["Name_score"][0], 1.0000)
      self.assertEqual(res["Name_score"][1], 1.0000)

    def test_1toManyMatchingOnDB(self):
      metaboliteA: MeMoMetabolite = MeMoMetabolite()
      metaboliteB: MeMoMetabolite = MeMoMetabolite()
      metaboliteA.set_annotations({"DatabaseA" : ["stuff","stuff3"]})
      metaboliteB.set_annotations({"DatabaseA" : ["stuff"]})
      metaboliteA.set_id("A")
      metaboliteB.set_id("B")
      model = MeMoModel([metaboliteA,metaboliteB])
      model2 = MeMoModel([metaboliteA, metaboliteB])

      res = model.match(model2, keepAllMatches = True)
      self.assertEqual(res.shape[0], 4)
      self.assertEqual(res["DB_score"][0], 1.0000)
      val = res["DB_score"][1]
      self.assertTrue(math.isclose(val, 0.5, rel_tol=1e-2))

      res = model.match(model2, keepAllMatches = False)
      self.assertEqual(res.shape[0], 2)
      self.assertEqual(res["DB_score"][0], 1.0000)
      self.assertEqual(res["DB_score"][1], 1.0000)

      res = model.match(model2, keepAllMatches = True, threshold_DB = 0.6)
      self.assertEqual(res.shape[0], 2)
      self.assertEqual(res["DB_score"][0], 1.0000)
      self.assertEqual(res["DB_score"][1], 1.0000)

    def test_1toManyMatchingOnInchi(self):
      metaboliteA: MeMoMetabolite = MeMoMetabolite()
      metaboliteB: MeMoMetabolite = MeMoMetabolite()
      metaboliteA.set_id("A")
      metaboliteB.set_id("B")
      metaboliteA.set_inchi_string("InChI=1S/H2O/h1H2")
      metaboliteB.set_inchi_string("InChI=1S/CH4/h1H4")

      model = MeMoModel([metaboliteA, metaboliteB])
      model2 = MeMoModel([metaboliteA, metaboliteB])
      res = model.match(model2, keepAllMatches = True)
      self.assertEqual(res.shape[0], 2)
      self.assertEqual(res["inchi_score"][0], 1.0000)
      val = res["inchi_score"][1]
      self.assertTrue(val==1)

      res2 = model.match(model2, keepAllMatches = False)
      self.assertTrue(all(res==res2)) # it does not make sense to have differences in the 1toMany cases for the inchis

    def test_keepUnmatched(self):
      metaboliteA: MeMoMetabolite = MeMoMetabolite()
      metaboliteB: MeMoMetabolite = MeMoMetabolite()
      metaboliteC: MeMoMetabolite = MeMoMetabolite()
      metaboliteA.set_id("A")
      metaboliteB.set_id("B")
      metaboliteC.set_id("C")
      metaboliteA.set_inchi_string("InChI=1S/H2O/h1H2")
      metaboliteB.set_inchi_string("InChI=1S/CH4/h1H4")
      metaboliteC.set_inchi_string(None)

      model = MeMoModel([metaboliteA, metaboliteC])
      model2 = MeMoModel([metaboliteA, metaboliteB])
      res = model.match(model2, keepAllMatches = False, keepUnmatched=True)
      self.assertEqual(res.shape[0], 3)
      self.assertEqual(pd.notna(res["met_id2"]).iloc[1], False)
      self.assertEqual(pd.notna(res["met_id1"]).iloc[2], False)

    def test_annotationCount(self):
      #this_directory = Path(__file__).parent
      #dat = this_directory.joinpath("../manually_merged_models")
      mod = cb.io.load_model("textbook")
      mod = MeMoModel.fromModel(mod)
      le = len(mod.metabolites)
      print(f"Amount of metabs {le}")
      print(f"Amount of unannotated inchis {sum([x._inchi_string == None for x in mod.metabolites])}")
      print(f"Annoatted {mod.annotate()}")
      print(f"Amount of unannotated inchis after Annotation {sum([x._inchi_string == None for x in mod.metabolites])}")

class Test_removeDuplicates(unittest.TestCase):
    this_directory = Path(__file__).parent
    dat = this_directory.joinpath("dat")
    
    def test_detectDuplicates(self):
        # load data
        mod_ori = cb.io.read_sbml_model(self.dat.joinpath("tiny_ecoli_keep_inchi.xml"))
        mmm_ori = MeMoModel.fromModel(mod_ori)
        mod_dup = cb.io.read_sbml_model(self.dat.joinpath("tiny_ecoli_keep_inchi_withduplicates.xml"))
        mmm_dup = MeMoModel.fromModel(mod_dup)
        mmm_dup_copy = mmm_dup.copy()
        # first make sure that the testcase model is different to the model with the duplicated 
        self.assertFalse(mmm_dup == mmm_ori)
        self.assertTrue(mmm_dup == mmm_dup_copy)

        # find duplicates
        dups = mmm_dup.match(mmm_dup)
        dups = detectDuplicates(dups)
        self.assertTrue(dups == [('pyr','pyruv')])

        # remove duplicates
        mmm_dup,rm_results = removeDuplicateMetabolites(mmm_dup)
        self.assertFalse(mmm_dup == mmm_dup_copy)
        self.assertTrue(mmm_dup == mmm_ori)

if __name__ == '__main__':
    unittest.main()
