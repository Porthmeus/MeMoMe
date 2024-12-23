import unittest
from pathlib import Path
import pandas as pd
import os
from unittest.mock import patch
from io import *
from src.annotation.annotateModelSEED import annotateModelSEED, annotateModelSEED_id, annotateModelSEED_entry
from src.annotation.annotateChEBI import annotateChEBI
from src.annotation.annotateBiGG import annotateBiGG, annotateBiGG_id, annotateBiGG_entry, handle_bigg_entries
from src.annotation.annotateVMH import annotateVMH_entry, annotateVMH, annotateVMH_id
from src.annotation.annotateHmdb import annotateHMDB_entry, annotateHMDB
from src.annotation.annotateModelSEED import annotateModelSEED_entry, extractModelSEEDAnnotationsFromAlias, extractModelSEEDAnnotationsFromAlias2
from src.annotation.annotateAux import AnnotationResult
from src.MeMoMetabolite import MeMoMetabolite


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
  def testBiggIDs(self, mock_err):
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


  def testBiggEntryToy(self):

    data = {
      'universal_bigg_id': ["1", "2", "3"],
      'name': ["Metabolite 1", "Metabolite 2", "Metabolite 3"],
      'database_links': ["http://identifiers.org/A/A_M1", "http://identifiers.org/A/A_M2" ,"http://identifiers.org/A/A_M3"]
    }
    
    # Create DataFrame
    df = pd.DataFrame(data)
    ret = annotateBiGG_entry("1", df, allow_missing_dbs = False)
    exprected = ({'A': ['A_M1']}, ['Metabolite 1'])
    self.assertEqual(ret, exprected)
  
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


  def testHMDBEntry(self):
    this_directory = Path(__file__).parent
    dbs_dir = this_directory.parent/Path("Databases")
    ret = annotateHMDB_entry("HMDB00972", allow_missing_dbs = False)
    self.assertTrue(len(ret[0]) > 0)
    self.assertTrue(len(ret[1]) > 0)
    
    ret = annotateHMDB_entry("", allow_missing_dbs = False)
    self.assertEqual(ret, (dict(), list()))


  def testSEEDEntry(self):
    self.maxDiff = None
    this_directory = Path(__file__).parent
    dbs_dir = this_directory.parent/Path("Databases")
    ret = annotateModelSEED_entry("cpd00052", allow_missing_dbs = False)
    self.assertEqual(sorted(ret[1]), sorted(["cytidine-triphosphate",
   "Cytidine triphosphate",
   "Cytidine 5'-triphosphate",
   "cytidine-5'-triphosphate",
   "CTP"]))
    print(ret[0])
    self.assertFalse(len(ret[0]) == 0)

    ret = annotateModelSEED_entry("", allow_missing_dbs = False)
    self.assertEqual(ret, (dict(), list()))



class Test_annotateID(unittest.TestCase):
  def testBiggID(self):
    this_directory = Path(__file__).parent
    dbs_dir = this_directory.parent/Path("Databases")
    metabolite: MeMoMetabolite = MeMoMetabolite()
    metabolite.set_id("13dpg")
    metabolite.set_names(["A"])
    
    ret = annotateBiGG_id([metabolite], allow_missing_dbs = False)
    expected_annotations = {'reactome': ['R-ALL-29800'], 'kegg.compound': ['C00236'], 'chebi': ['CHEBI:11881', 'CHEBI:16001', 'CHEBI:1658', 'CHEBI:20189', 'CHEBI:57604'], 'hmdb': ['HMDB62758'], 'inchikey': ['LJQLQCAXBUHEAZ-UWTATZPHSA-J'], 'biocyc': ['META:DPG'], 'metanetx.chemical': ['MNXM261'], 'seed.compound': ['cpd00203']}
    self.assertEqual(metabolite.annotations, expected_annotations)
    self.assertEqual(metabolite.names, ["3-Phospho-D-glyceroyl phosphate", "A"])
    self.assertEqual(ret, AnnotationResult(0, 1, 1))


  def testVMH_id(self):
    this_directory = Path(__file__).parent
    dbs_dir = this_directory.parent/Path("Databases")
    metabolite: MeMoMetabolite = MeMoMetabolite()
    metabolite.set_id("10fthf")
    metabolite.set_names(["A"])
    ret = annotateVMH_id([metabolite], allow_missing_dbs = False)

    expected_annotations =  {'biggId': ['10fthf'], 'keggId': ['C00234'], 'cheBlId': ['15637'], 'inchiString': ['InChI=1S/C20H23N7O7/c21-20-25-16-15(18(32)26-20)23-11(7-22-16)8-27(9-28)12-3-1-10(2-4-12)17(31)24-13(19(33)34)5-6-14(29)30/h1-4,9,11,13,23H,5-8H2,(H,24,31)(H,29,30)(H,33,34)(H4,21,22,25,26,32)/p-2/t11?,13-/m0/s1'], 'inchiKey': ['AUFGTPPARQZWDO-YUZLPWPTSA-L'], 'smile': ['[H]OC1=NC(=NC2=C1N([H])C([H])(C([H])([H])N(C([H])=O)C1=C([H])C([H])=C(C([H])=C1[H])C(=O)N([H])[C@]([H])(C([O-])=O)C([H])([H])C([H])([H])C([O-])=O)C([H])([H])N2[H])N([H])[H]'], 'hmdb': ['HMDB0000972'], 'metanetx': ['MNXM237'], 'seed': ['cpd00201'], 'biocyc': ['10-FORMYL-THF']}
    self.assertEqual(metabolite.annotations, expected_annotations)
    self.assertEqual(metabolite.names, ['10-Formyltetrahydrofolate', 'A'])
    self.assertEqual(ret, AnnotationResult(0, 1, 1))


class Test_annotateFull(unittest.TestCase):
  def testBiggAnnotate(self):
    metabolite: MeMoMetabolite = MeMoMetabolite()
    metabolite.set_id("13dpg")
    metabolite.annotations = {'bigg.metabolite': ['13dpg']}
    
    ret = annotateBiGG([metabolite], allow_missing_dbs = False)
    expected_annotations = {'bigg.metabolite': ['13dpg'], 'reactome': ['R-ALL-29800'], 'kegg.compound': ['C00236'], 'chebi': ['CHEBI:11881', 'CHEBI:16001', 'CHEBI:1658', 'CHEBI:20189', 'CHEBI:57604'], 'hmdb': ['HMDB62758'], 'inchikey': ['LJQLQCAXBUHEAZ-UWTATZPHSA-J'], 'biocyc': ['META:DPG'], 'metanetx.chemical': ['MNXM261'], 'seed.compound': ['cpd00203']}
    self.assertEqual(metabolite.annotations, expected_annotations)
    self.assertEqual(metabolite.names, ["3-Phospho-D-glyceroyl phosphate"])
    self.assertEqual(ret, AnnotationResult(0, 1, 1))


  def testVMHAnnotate(self):
    metabolite: MeMoMetabolite = MeMoMetabolite()
    metabolite.set_id("10fthf")
    metabolite.annotations = {'vmhmetabolite': ['10fthf']}
    ret = annotateVMH([metabolite], allow_missing_dbs = False)

    expected_annotations =  {'vmhmetabolite': ['10fthf'], 'biggId': ['10fthf'], 'keggId': ['C00234'], 'cheBlId': ['15637'], 'inchiString': ['InChI=1S/C20H23N7O7/c21-20-25-16-15(18(32)26-20)23-11(7-22-16)8-27(9-28)12-3-1-10(2-4-12)17(31)24-13(19(33)34)5-6-14(29)30/h1-4,9,11,13,23H,5-8H2,(H,24,31)(H,29,30)(H,33,34)(H4,21,22,25,26,32)/p-2/t11?,13-/m0/s1'], 'inchiKey': ['AUFGTPPARQZWDO-YUZLPWPTSA-L'], 'smile': ['[H]OC1=NC(=NC2=C1N([H])C([H])(C([H])([H])N(C([H])=O)C1=C([H])C([H])=C(C([H])=C1[H])C(=O)N([H])[C@]([H])(C([O-])=O)C([H])([H])C([H])([H])C([O-])=O)C([H])([H])N2[H])N([H])[H]'], 'hmdb': ['HMDB0000972'], 'metanetx': ['MNXM237'], 'seed': ['cpd00201'], 'biocyc': ['10-FORMYL-THF']}
    self.assertEqual(metabolite.annotations, expected_annotations)
    self.assertEqual(metabolite.names, ['10-Formyltetrahydrofolate'])
    self.assertEqual(ret, AnnotationResult(0, 1, 1))

  def testChEBIAnnotate(self):
    metabolite: MeMoMetabolite = MeMoMetabolite()
    metabolite.set_id('117228')
    metabolite.annotations = {'chebi': ['117228']}
    ret = annotateChEBI([metabolite], allow_missing_dbs = False)
    expected_inchi = "InChI=1S/C20H27N3O2/c1-4-20(3)13-14-9-6-7-10-15(14)17-16(20)18(25)23(5-2)19(22-17)21-11-8-12-24/h6-7,9-10,24H,4-5,8,11-13H2,1-3H3,(H,21,22)"
    self.assertEqual(ret, AnnotationResult(1, 0, 0))
    self.assertEqual(metabolite._inchi_string, expected_inchi)

  def testHMDBAnnotate(self):
    metabolite: MeMoMetabolite = MeMoMetabolite()
    metabolite.set_id("10fthf")
    metabolite.annotations = {'HMDB': ['HMDB0000972']}
    expected_annotations = {'HMDB': ['HMDB0000972'], 'accession': ['HMDB0000972'], 'chemical_formula': ['C20H23N7O7'], 'iupac_name': ['(2S)-2-[(4-{N-[(4-hydroxy-2-imino-1,2,5,6,7,8-hexahydropteridin-6-yl)methyl]formamido}phenyl)formamido]pentanedioic acid'], 'traditional_iupac': ['(2S)-2-[(4-{N-[(4-hydroxy-2-imino-5,6,7,8-tetrahydro-1H-pteridin-6-yl)methyl]formamido}phenyl)formamido]pentanedioic acid'], 'smiles': ['NC1=NC(=O)C2=C(NCC(CN(C=O)C3=CC=C(C=C3)C(=O)N[C@@H](CCC(O)=O)C(O)=O)N2)N1'], 'inchi': ['InChI=1S/C20H23N7O7/c21-20-25-16-15(18(32)26-20)23-11(7-22-16)8-27(9-28)12-3-1-10(2-4-12)17(31)24-13(19(33)34)5-6-14(29)30/h1-4,9,11,13,23H,5-8H2,(H,24,31)(H,29,30)(H,33,34)(H4,21,22,25,26,32)/t11?,13-/m0/s1'], 'chebi_id': [15637.0], 'pubchem_compound_id': [122347.0], 'kegg_id': ['C00234'], 'bigg_id': [34337.0], 'vmh_id': ['10FTHF']}
    ret = annotateHMDB([metabolite], allow_missing_dbs = False)
    self.assertEqual(metabolite.annotations, expected_annotations)
    self.assertEqual(metabolite.names, ['10-Formyltetrahydrofolate'])
    self.assertEqual(ret, AnnotationResult(0, 1, 1))


class Test_annotateAuxiliares(unittest.TestCase):

  this_directory = Path(__file__).parent
  dbs_dir = this_directory.parent/Path("Databases")

  def test_extractModelSEEDAnnotationsFromAlias(self):
    # Aliases of cpd00001
    aliases1 = "Name: H20; H2O; H3O+; HO-; Hydroxide ion; OH; OH-; Water; hydrogen oxide; hydroxide; hydroxide ion; hydroxyl; hydroxyl ion; oxonium; water|AraCyc: OH; WATER|BiGG: h2o; oh1|BrachyCyc: WATER|KEGG: C00001; C01328|MetaCyc: OH; OXONIUM; WATER"

    aliases2 = "Name: NADP(H); NADP-red; NADP-reduced; NADPH; NADPH+H+; NADPH2; Nicotinamide adenine dinucleotide phosphate - reduced; Nicotinamide adenine dinucleotide phosphate-reduced; Nicotinamideadeninedinucleotidephosphate-reduced; Reduced nicotinamide adenine dinucleotide phosphate; TPNH; beta-NADPH; dihydronicotinamide adenine dinucleotide phosphate; dihydronicotinamide adenine dinucleotide phosphate reduced; dihydronicotinamide adenine dinucleotide-P; dihydrotriphosphopyridine nucleotide; dihydrotriphosphopyridine nucleotide reduced; reduced NADP; reduced dihydrotriphosphopyridine nucleotide; reduced nicotinamide adenine dinucleotide phosphate|AraCyc: NADPH|BiGG: nadph|BrachyCyc: NADPH|KEGG: C00005|MetaCyc: NADPH"

    # Aliases of cpd00002
    aliases3 = "Name: ATP; Adenosine 5'-triphosphate; adenosine-5'-triphosphate; adenosine-triphosphate; adenylpyrophosphate|AraCyc: ATP|BiGG: atp|BrachyCyc: ATP|KEGG: C00002|MetaCyc: ATP"

    expected1 = ({'AraCyc': ['OH', 'WATER'], 'BiGG': ['h2o', 'oh1'], 'BrachyCyc': ['WATER'], 'KEGG': ['C00001', 'C01328'], 'MetaCyc': ['OH', 'OXONIUM', 'WATER']}, ['H20', 'H2O', 'H3O+', 'HO-', 'Hydroxide ion', 'OH', 'OH-', 'Water', 'hydrogen oxide', 'hydroxide', 'hydroxide ion', 'hydroxyl', 'hydroxyl ion', 'oxonium', 'water'])
    extracted1 = extractModelSEEDAnnotationsFromAlias(aliases1)
    extracted11 = extractModelSEEDAnnotationsFromAlias2(aliases1)
    print(extracted1)
    print(extracted11)
    #self.assertEqual(extracted1, expected1)

    #extracted35 = extractModelSEEDAnnotationsFromAlias2(aliases3)
    #expected3 = ({'AraCyc': ['ATP'], 'BiGG': ['atp'], 'BrachyCyc': ['ATP'], 'KEGG': ['C00002'], 'MetaCyc': ['ATP']}, ['ATP', "Adenosine 5'-triphosphate", "adenosine-5'-triphosphate", 'adenosine-triphosphate', 'adenylpyrophosphate'])
    #print(extracted3)
    #print(extracted35)
    ##self.assertEqual(extracted3, expected3)
