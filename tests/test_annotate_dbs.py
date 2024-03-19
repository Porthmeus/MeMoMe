import unittest
from pathlib import Path
import cobra as cb
import pandas as pd
import os
from unittest.mock import patch
from io import *
import warnings
from src.MeMoMetabolite import MeMoMetabolite
from src.MeMoModel import MeMoModel
from src.annotateModelSEED import annotateModelSEED, annotateModelSEED_id, correctAnnotationKeys, annotateModelSEED_entry
from src.annotateChEBI import annotateChEBI
from src.annotateBiGG import annotateBiGG, annotateBiGG_id, annotateBiGG_entry
from src.annotateAux import AnnotationResult


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

    self.makeDbVis("BiGG.tsvxxx")

  @patch('sys.stderr', new_callable=StringIO)
  def testChEBI(self, mock_err):
    self.makeDbInvis("chebiId_inchi.tsv")
    self.assertEqual(annotateChEBI([], allow_missing_dbs = True), AnnotationResult(0,0,0))
    output = mock_err.getvalue().strip()
    self.assertIn("No such file or directory", output)
    
    with self.assertRaises(FileNotFoundError):
      annotateChEBI([],  allow_missing_dbs = False)

    self.makeDbVis("chebiId_inchi.tsvxxx")

  @patch('sys.stderr', new_callable=StringIO)
  def testSEED(self, mock_err):
    self.makeDbInvis("modelSeed.tsv")
    self.assertEqual(annotateModelSEED([], allow_missing_dbs = True), AnnotationResult(0,0,0))
    output = mock_err.getvalue().strip()
    self.assertIn("No such file or directory", output)
    
    with self.assertRaises(FileNotFoundError):
      annotateModelSEED([],  allow_missing_dbs = False)

    self.makeDbVis("modelSeed.tsvxxx")

  @patch('sys.stderr', new_callable=StringIO)
  def testSEED_id(self, mock_err):
    self.makeDbInvis("modelSeed.tsv")
    self.assertEqual(annotateModelSEED_id([], allow_missing_dbs = True), AnnotationResult(0,0,0))
    output = mock_err.getvalue().strip()
    self.assertIn("No such file or directory", output)
    
    with self.assertRaises(FileNotFoundError):
      annotateModelSEED_id([],  allow_missing_dbs = False)

    self.makeDbVis("modelSeed.tsvxxx")

  @patch('sys.stderr', new_callable=StringIO)
  def testSEED_entry(self, mock_err):
    self.makeDbInvis("modelSeed.tsv")
    self.assertEqual(annotateModelSEED_entry("", pd.DataFrame(), allow_missing_dbs = True), ({}, [], {}, {}))
    output = mock_err.getvalue().strip()
    self.assertIn("No such file or directory", output)
    
    with self.assertRaises(FileNotFoundError):
      annotateModelSEED_entry("", pd.DataFrame(),  allow_missing_dbs = False)

    self.makeDbVis("modelSeed.tsvxxx")
  @patch('sys.stderr', new_callable=StringIO)
  def correctAnnotationKeys(self, mock_err):
    self.makeDbInvis("modelSeed.tsv")
    self.assertEqual(correctAnnotationKeys({}, allow_missing_dbs = True), AnnotationResult(0,0,0))
    output = mock_err.getvalue().strip()
    self.assertIn("No such file or directory", output)
    
    with self.assertRaises(FileNotFoundError):
      correctAnnotationKeys({},  allow_missing_dbs = False)

    self.makeDbVis("modelSeed.tsvxxx")
