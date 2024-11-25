import unittest
from pathlib import Path
import pandas as pd
import os
from unittest.mock import patch
from io import *
from src.annotation.annotateModelSEED import annotateModelSEED, annotateModelSEED_id, correctAnnotationKeys, annotateModelSEED_entry
from src.annotation.annotateChEBI import annotateChEBI
from src.annotation.annotateBiGG import annotateBiGG, annotateBiGG_id, annotateBiGG_entry, handle_bigg_entries
from src.annotation.annotateVMH import annotateVMH_entry, annotateVMH, annotateVMH_id
from src.annotation.annotateAux import AnnotationResult


# List of files to check for
files_to_check = ["BiGG.tsv", "chebiId_inchi.tsv", "identifiers.org_registry.json", "modelSeed.tsv", "vmh.json"]



# Check if all files are present in the directory

class Test_annotateMissingDbs(unittest.TestCase):
  # The directory of this file
  #this_directory = Path("tests")
  this_directory = Path(__file__).parent
  dbs_dir = this_directory.parent/Path("Databases")



  def makeDbInvis(self, old_name):
    """ 
    This functions adds a xxx to the database name to make to temporarily not available
    old_name: Name of a database without xxx as a suffix
    return: Name of a database that ends in xxx
    """
    # Construct the new filename with "xxx" appended
    new_name = f"{old_name}xxx"
    old_path = os.path.join(self.dbs_dir, old_name)
    new_path = os.path.join(self.dbs_dir, new_name)
    os.rename(old_path, new_path)
  
  def makeDbVis(self, new_name):
    """ 
    This functions removes a xxx to the database name to make it available again
    new_name: Name of a database that ends in xxx
    return: Name of a database without xxx as a suffix
    """
    original_name = new_name.rsplit("xxx", 1)[0] 
      
    old_path = os.path.join(self.dbs_dir, original_name)
    new_path = os.path.join(self.dbs_dir, new_name)
    try: 
      os.rename(new_path, old_path) 
    except Exception as e:
      pass

  def setUp(self):
    # Check if all databases are there before we do a test
    self.assertTrue(all((self.dbs_dir/ file).is_file() for file in files_to_check))
  
  def tearDown(self):
    self.makeDbVis("BiGG.tsvxxx")
    self.makeDbVis("chebiId_inchi.tsvxxx")
    self.makeDbVis("modelSeed.tsvxxx")
    self.makeDbVis("vmh.jsonxxx")
  
  
  @patch('sys.stderr', new_callable=StringIO)
  def testBiggEntry(self, mock_err):
    self.makeDbInvis("BiGG.tsv")
    self.assertEqual(annotateBiGG_entry("", pd.DataFrame(), allow_missing_dbs = True), ({}, []))
    output = mock_err.getvalue().strip()
    self.assertIn("No such file or directory", output)
    
    with self.assertRaises(FileNotFoundError):
      annotateBiGG_entry("", pd.DataFrame(), allow_missing_dbs = False)

    self.makeDbVis("BiGG.tsvxxx")

  @patch('sys.stderr', new_callable=StringIO)
  def testBiggID(self, mock_err):
    self.makeDbInvis("BiGG.tsv")
    self.assertEqual(annotateBiGG_id([], allow_missing_dbs = True), AnnotationResult(0,0,0))
    output = mock_err.getvalue().strip()
    self.assertIn("No such file or directory", output)
    
    with self.assertRaises(FileNotFoundError):
      annotateBiGG_id([], allow_missing_dbs = False)

    self.makeDbVis("BiGG.tsvxxx")


  @patch('sys.stderr', new_callable=StringIO)
  def testBigg(self, mock_err):
    self.makeDbInvis("BiGG.tsv")
    self.assertEqual(annotateBiGG([], allow_missing_dbs = True), AnnotationResult(0,0,0))
    output = mock_err.getvalue().strip()
    self.assertIn("No such file or directory", output)
    
    with self.assertRaises(FileNotFoundError):
      annotateBiGG([], allow_missing_dbs = False)
      self.assertIn("No such file or directory", output)

    self.makeDbVis("BiGG.tsvxxx")

  @patch('sys.stderr', new_callable=StringIO)
  def testChEBI(self, mock_err):
    self.makeDbInvis("chebiId_inchi.tsv")
    self.assertEqual(annotateChEBI([], allow_missing_dbs = True), AnnotationResult(0,0,0))
    output = mock_err.getvalue().strip()
    self.assertIn("No such file or directory", output)
    
    with self.assertRaises(FileNotFoundError):
      annotateChEBI([],  allow_missing_dbs = False)
      self.assertIn("No such file or directory", output)

    self.makeDbVis("chebiId_inchi.tsvxxx")

  @patch('sys.stderr', new_callable=StringIO)
  def testSEED(self, mock_err):
    self.makeDbInvis("modelSeed.tsv")
    self.assertEqual(annotateModelSEED([], allow_missing_dbs = True), AnnotationResult(0,0,0))
    output = mock_err.getvalue().strip()
    self.assertIn("No such file or directory", output)
    
    with self.assertRaises(FileNotFoundError):
      annotateModelSEED([],  allow_missing_dbs = False)
      self.assertIn("No such file or directory", output)

    self.makeDbVis("modelSeed.tsvxxx")

  @patch('sys.stderr', new_callable=StringIO)
  def testSEED_id(self, mock_err):
    self.makeDbInvis("modelSeed.tsv")
    self.assertEqual(annotateModelSEED_id([], allow_missing_dbs = True), AnnotationResult(0,0,0))
    output = mock_err.getvalue().strip()
    self.assertIn("No such file or directory", output)
    
    with self.assertRaises(FileNotFoundError):
      annotateModelSEED_id([],  allow_missing_dbs = False)
      self.assertIn("No such file or directory", output)

    self.makeDbVis("modelSeed.tsvxxx")

  @patch('sys.stderr', new_callable=StringIO)
  def testSEED_entry(self, mock_err):
    self.makeDbInvis("modelSeed.tsv")
    self.assertEqual(annotateModelSEED_entry("", pd.DataFrame(), allow_missing_dbs = True), ({}, [], {}, {}))
    output = mock_err.getvalue().strip()
    self.assertIn("No such file or directory", output)
    
    with self.assertRaises(FileNotFoundError):
      annotateModelSEED_entry("", pd.DataFrame(),  allow_missing_dbs = False)
      self.assertIn("No such file or directory", output)

    self.makeDbVis("modelSeed.tsvxxx")
  
  @patch('sys.stderr', new_callable=StringIO)
  def correctAnnotationKeys(self, mock_err):
    self.makeDbInvis("modelSeed.tsv")
    self.assertEqual(correctAnnotationKeys({}, allow_missing_dbs = True), AnnotationResult(0,0,0))
    output = mock_err.getvalue().strip()
    self.assertIn("No such file or directory", output)
    
    with self.assertRaises(FileNotFoundError):
      correctAnnotationKeys({},  allow_missing_dbs = False)
      self.assertIn("No such file or directory", output)

    self.makeDbVis("modelSeed.tsvxxx")


  @patch('sys.stderr', new_callable=StringIO)
  def testVMH_entry(self, mock_err):
    self.makeDbInvis("vmh.json")
    self.assertEqual(annotateVMH_entry("", pd.DataFrame(), allow_missing_dbs = True), ({}, []))
    output = mock_err.getvalue().strip()
    self.assertIn("No such file or directory", output)
    
    with self.assertRaises(FileNotFoundError):
      annotateVMH_entry("", pd.DataFrame(),  allow_missing_dbs = False)
      self.assertIn("No such file or directory", output)
    self.makeDbVis("vmh.json")



  @patch('sys.stderr', new_callable=StringIO)
  def testVMH(self, mock_err):
    self.makeDbInvis("vmh.json")
    self.assertEqual(annotateVMH([], allow_missing_dbs = True),AnnotationResult(0, 0, 0))
    output = mock_err.getvalue().strip()
    self.assertIn("No such file or directory", output)
    
    with self.assertRaises(FileNotFoundError):
      annotateVMH([],  allow_missing_dbs = False)
      self.assertIn("No such file or directory", output)
    self.makeDbVis("vmh.json")

  @patch('sys.stderr', new_callable=StringIO)
  def testVMH_id(self, mock_err):
    self.makeDbInvis("vmh.json")
    self.assertEqual(annotateVMH_id([], allow_missing_dbs = True),AnnotationResult(0, 0, 0))
    output = mock_err.getvalue().strip()
    self.assertIn("No such file or directory", output)
    
    with self.assertRaises(FileNotFoundError):
      annotateVMH_id([],  allow_missing_dbs = False)
      self.assertIn("No such file or directory", output)
    self.makeDbVis("vmh.json")



class Test_annotateEntryFunctions(unittest.TestCase):

  def testBiggEntry(self):
    this_directory = Path(__file__).parent
    dbs_dir = this_directory.parent/Path("Databases")
    ret = annotateBiGG_entry("13dpg", allow_missing_dbs = False)
    self.assertTrue(len(ret[0]) > 0)
    self.assertTrue(len(ret[1]) > 0)
    
    ret = annotateBiGG_entry("", allow_missing_dbs = False)
    self.assertEqual(ret, (dict(), list()))
  
  def testBiggHandler(self):
    urls = pd.Series([
      "http://identifiers.org/hmdb/HMDB02322",
      "http://identifiers.org/hmdb/HMDB04567",
      "http://identifiers.org/hmdb/HMDB07890",
      "http://identifiers.org/A/A_M1"
    ])

    res = handle_bigg_entries(urls)
    self.assertTrue("hmdb" in res)
    self.assertTrue("A" in res)

    self.assertTrue(set(["HMDB02322", "HMDB04567", "HMDB07890"]) == set(res["hmdb"]))
    self.assertTrue(set(["A_M1"]) == set(res["A"]))


 
  def testVMHEntry(self):
    this_directory = Path(__file__).parent
    dbs_dir = this_directory.parent/Path("Databases")
    ret = annotateVMH_entry("10fthf", allow_missing_dbs = False)
    
    self.assertEqual(ret[1], ["10-Formyltetrahydrofolate"])
    self.assertFalse(len(ret[0]) == 0)
