import unittest
from pathlib import Path
import pandas as pd
import os
from unittest.mock import patch
from io import *
from src.annotation.annotateAux import AnnotationResult, AnnotationKey, annotateEntry, load_database, handleIDs, handleMetabolites
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
  
  
  #@patch('sys.stderr', new_callable=StringIO)
  #def testBiggEntry(self, mock_err):
  #  self.makeDbInvis("BiGG.tsv")
  #  self.assertEqual(annotateBiGG_entry("", pd.DataFrame(), allow_missing_dbs = True), ({}, []))
  #  output = mock_err.getvalue().strip()
  #  self.assertIn("No such file or directory", output)
  #  with self.assertRaises(FileNotFoundError):
  #    annotateBiGG_entry("", pd.DataFrame(), allow_missing_dbs = False)

  #  self.makeDbVis("BiGG.tsvxxx")

  #@patch('sys.stderr', new_callable=StringIO)
  #def testBiggIDs(self, mock_err):
  #  self.makeDbInvis("BiGG.tsv")
  #  self.assertEqual(annotateBiGG_id([], allow_missing_dbs = True), AnnotationResult(0,0,0))
  #  output = mock_err.getvalue().strip()
  #  self.assertIn("No such file or directory", output)
  #  
  #  with self.assertRaises(FileNotFoundError):
  #    annotateBiGG_id([], allow_missing_dbs = False)

  #  self.makeDbVis("BiGG.tsvxxx")


  #@patch('sys.stderr', new_callable=StringIO)
  #def testBigg(self, mock_err):
  #  self.makeDbInvis("BiGG.tsv")
  #  self.assertEqual(annotateBiGG([], allow_missing_dbs = True), AnnotationResult(0,0,0))
  #  output = mock_err.getvalue().strip()
  #  self.assertIn("No such file or directory", output)
  #  
  #  with self.assertRaises(FileNotFoundError):
  #    annotateBiGG([], allow_missing_dbs = False)
  #    self.assertIn("No such file or directory", output)

  #  self.makeDbVis("BiGG.tsvxxx")

  #@patch('sys.stderr', new_callable=StringIO)
  #def testChEBI(self, mock_err):
  #  self.makeDbInvis("chebiId_inchi.tsv")
  #  self.assertEqual(annotateChEBI([], allow_missing_dbs = True), AnnotationResult(0,0,0))
  #  output = mock_err.getvalue().strip()
  #  self.assertIn("No such file or directory", output)
  #  
  #  with self.assertRaises(FileNotFoundError):
  #    annotateChEBI([],  allow_missing_dbs = False)
  #    self.assertIn("No such file or directory", output)

  #  self.makeDbVis("chebiId_inchi.tsvxxx")

  #@patch('sys.stderr', new_callable=StringIO)
  #def testSEED(self, mock_err):
  #  self.makeDbInvis("modelSeed.tsv")
  #  self.assertEqual(annotateModelSEED([], allow_missing_dbs = True), AnnotationResult(0,0,0))
  #  output = mock_err.getvalue().strip()
  #  self.assertIn("No such file or directory", output)
  #  
  #  with self.assertRaises(FileNotFoundError):
  #    annotateModelSEED([],  allow_missing_dbs = False)
  #    self.assertIn("No such file or directory", output)

  #  self.makeDbVis("modelSeed.tsvxxx")

  #@patch('sys.stderr', new_callable=StringIO)
  #def testSEED_id(self, mock_err):
  #  self.makeDbInvis("modelSeed.tsv")
  #  self.assertEqual(annotateModelSEED_id([], allow_missing_dbs = True), AnnotationResult(0,0,0))
  #  output = mock_err.getvalue().strip()
  #  self.assertIn("No such file or directory", output)
  #  
  #  with self.assertRaises(FileNotFoundError):
  #    annotateModelSEED_id([],  allow_missing_dbs = False)
  #    self.assertIn("No such file or directory", output)

  #  self.makeDbVis("modelSeed.tsvxxx")

  #@patch('sys.stderr', new_callable=StringIO)
  #def testSEED_entry(self, mock_err):
  #  self.makeDbInvis("modelSeed.tsv")
  #  self.assertEqual(annotateModelSEED_entry("", pd.DataFrame(), allow_missing_dbs = True), ({}, []))
  #  output = mock_err.getvalue().strip()
  #  self.assertIn("No such file or directory", output)
  #  
  #  with self.assertRaises(FileNotFoundError):
  #    annotateModelSEED_entry("", pd.DataFrame(),  allow_missing_dbs = False)
  #    self.assertIn("No such file or directory", output)

  #  self.makeDbVis("modelSeed.tsvxxx")
  
#  @patch('sys.stderr', new_callable=StringIO)
#  def correctAnnotationKeys(self, mock_err):
#    self.makeDbInvis("modelSeed.tsv")
#    self.assertEqual(correctAnnotationKeys({}, allow_missing_dbs = True), AnnotationResult(0,0,0))
#    output = mock_err.getvalue().strip()
#    self.assertIn("No such file or directory", output)
#    
#    with self.assertRaises(FileNotFoundError):
#      correctAnnotationKeys({},  allow_missing_dbs = False)
#      self.assertIn("No such file or directory", output)
#
#    self.makeDbVis("modelSeed.tsvxxx")


  #@patch('sys.stderr', new_callable=StringIO)
  #def testVMH_entry(self, mock_err):
  #  self.makeDbInvis("vmh.json")
  #  self.assertEqual(annotateVMH_entry(AnnotationKey(""), pd.DataFrame(), allow_missing_dbs = True), ({}, []))
  #  output = mock_err.getvalue().strip()
  #  self.assertIn("No such file or directory", output)
  #  
  #  with self.assertRaises(FileNotFoundError):
  #    annotateVMH_entry(AnnotationKey(""), pd.DataFrame(),  allow_missing_dbs = False)
  #    self.assertIn("No such file or directory", output)
  #  self.makeDbVis("vmh.json")



  #@patch('sys.stderr', new_callable=StringIO)
  #def testVMH(self, mock_err):
  #  self.makeDbInvis("vmh.json")
  #  self.assertEqual(annotateVMH([], allow_missing_dbs = True),AnnotationResult(0, 0, 0))
  #  output = mock_err.getvalue().strip()
  #  self.assertIn("No such file or directory", output)
  #  
  #  with self.assertRaises(FileNotFoundError):
  #    annotateVMH([],  allow_missing_dbs = False)
  #    self.assertIn("No such file or directory", output)
  #  self.makeDbVis("vmh.json")

  #@patch('sys.stderr', new_callable=StringIO)
  #def testVMH_id(self, mock_err):
  #  self.makeDbInvis("vmh.json")
  #  self.assertEqual(annotateVMH_id([], allow_missing_dbs = True),AnnotationResult(0, 0, 0))
  #  output = mock_err.getvalue().strip()
  #  self.assertIn("No such file or directory", output)
  #  
  #  with self.assertRaises(FileNotFoundError):
  #    annotateVMH_id([],  allow_missing_dbs = False)
  #    self.assertIn("No such file or directory", output)
  #  self.makeDbVis("vmh.json")



class Test_annotateEntryFunctions(unittest.TestCase):

  def testBiggEntry(self):
    #ret = annotateBiGG_entry("13dpg", allow_missing_dbs = False)
    db = load_database("BiGG")
    ret = annotateEntry("13dpg", db)
    print(ret)
    self.assertTrue(len(ret[0]) > 0)
    self.assertTrue(len(ret[1]) > 0)
    
    ret = annotateEntry("", db)
    self.assertEqual(ret, (dict(), list(), ""))


  def testBiggEntryToy(self):

    data = {
      'id': ["1", "2", "3"],
      'name': ["Metabolite 1", "Metabolite 2", "Metabolite 3"],
      'DBs': ['''{"DB_A":["A_M1"]}''',
              '''{"DB_C":["A_M2", "A_M2"], "DB_D":["A_M3"]}''',
              '''{"DB_F":["A_M1", "A_M2"], "DB_G":["A_M3"]}'''],
      'inchi': ["InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1", "InChI=1S/C3H8O10P2/c4-2(1-12-14(6,7)8)3(5)13-15(9,10)11/h2,4H,1H2,(H2,6,7,8)(H2,9,10,11)/p-4/t2-/m1/s1","InChI=1S/C3H8O10P2/c4-2(1-12-14(6,7)8)3(5)13-1"]
    }
    
    # Create DataFrame
    df = pd.DataFrame(data)
    ret = annotateEntry("1", df)
    expected = ({'DB_A': ['A_M1']}, ['Metabolite 1'], "InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1")
    self.assertEqual(ret, expected)
  
#  def testBiggHandler(self):
#    urls = pd.Series([
#      "http://identifiers.org/hmdb/HMDB02322",
#      "http://identifiers.org/hmdb/HMDB04567",
#      "http://identifiers.org/hmdb/HMDB07890",
#      "http://identifiers.org/A/A_M1"
#    ])
#
#    res = handle_bigg_entries(urls)
#    self.assertTrue("hmdb" in res)
#    self.assertTrue("A" in res)
#
#    self.assertTrue(set(["HMDB02322", "HMDB04567", "HMDB07890"]) == set(res["hmdb"]))
#    self.assertTrue(set(["A_M1"]) == set(res["A"]))


 
  def testVMHEntry(self):
    db = load_database("VMH")
    ret = annotateEntry("10fthf", db)
    print(ret)
    self.assertEqual(ret[1][0], "10-Formyltetrahydrofolate")
    self.assertFalse(len(ret[0]) == 0)


  def testHMDBEntry(self):
    db = load_database("HMDB")
    ret = annotateEntry("HMDB00972", db)
    print(ret)
    self.assertTrue(len(ret[0]) > 0)
    self.assertTrue(len(ret[1]) > 0)
    
    ret = annotateEntry("", db)
    self.assertEqual(ret, (dict(), list(),""))

  def testSEEDEntry(self):
    db = load_database("ModelSeed")
    ret = annotateEntry("cpd00052", db)
    self.assertEqual(sorted(ret[1]), sorted(["cytidine-triphosphate",
   "cytidine triphosphate",
   "cytidine 5'-triphosphate",
   "cytidine-5'-triphosphate",
   "ctp"]))
    self.assertFalse(len(ret[0]) == 0)

    ret = annotateEntry("", db)
    self.assertEqual(ret, (dict(), list(), ""))



class Test_annotateID(unittest.TestCase):
  def testBiggID(self):
    this_directory = Path(__file__).parent
    dbs_dir = this_directory.parent/Path("Databases")
    metabolite: MeMoMetabolite = MeMoMetabolite()
    metabolite.set_id("13dpg")
    metabolite.set_names(["A"], source = "test")
    
    ret = handleIDs([metabolite], db_name = "BiGG", allow_missing_dbs = False)
    expected_annotations = {'reactome': ['R-ALL-29800'], 'kegg.compound': ['C00236'], 'chebi': ['CHEBI:11881', 'CHEBI:16001', 'CHEBI:1658', 'CHEBI:20189', 'CHEBI:57604'], 'hmdb': ['HMDB62758'], 'inchikey': ['LJQLQCAXBUHEAZ-UWTATZPHSA-J'], 'biocyc': ['META:DPG'], 'metanetx.chemical': ['MNXM261'], 'seed.compound': ['cpd00203']}
    self.assertEqual(metabolite.annotations, expected_annotations)
    self.assertEqual(metabolite.names, ["3-Phospho-D-glyceroyl phosphate", "A"])
    self.assertEqual(ret, AnnotationResult(1, 1, 1))


  def testVMH_id(self):
    this_directory = Path(__file__).parent
    dbs_dir = this_directory.parent/Path("Databases")
    metabolite: MeMoMetabolite = MeMoMetabolite()
    metabolite.set_id("10fthf")
    metabolite.set_names(["A"], source = "test")
    ret = handleIDs([metabolite], db_name = "VMH", allow_missing_dbs = False)
    self.maxDiff = None
    expected_annotations =  {'vmhmetabolite': ['10fthf'], 'bigg.metabolite': ['10fthf'], 'kegg.compound': ['C00234'], 'chemspider': ['109092'], 'chebi': ['15637'], 'biocyc': ['10-FORMYL-THF'], 'foodb.compound': ['FDB022345'], 'hmdb': ['HMDB0000972'], 'metanetx': ['MNXM237'], 'metlin': ['5912'], 'pubchem.compound': ['122347'], 'seed.compound': ['cpd00201'], 'cas': ['2800-34-2']}
    self.assertEqual(metabolite.annotations, expected_annotations)
    self.assertEqual(metabolite.names, ['(2S)-2-[(4-{N-[(2-amino-4-oxo-1,4,5,6,7,8-hexahydropteridin-6-yl)methyl]formamido}phenyl)formamido]pentanedioic acid','10-Formyltetrahydrofolate', 'A'])
    self.assertEqual(ret, AnnotationResult(1, 1, 1))


  def testSEED_id(self):
    this_directory = Path(__file__).parent
    dbs_dir = this_directory.parent/Path("Databases")
    metabolite: MeMoMetabolite = MeMoMetabolite()
    metabolite.set_id("cpd00002")
    ret = handleIDs([metabolite],db_name = "ModelSeed", allow_missing_dbs = False)
    expected_annotations = {'bigg.metabolite': ['atp'], 'kegg.compound': ['C00002'], 'metacyc.compound': ['ATP']}
    expected_names = ["adenosine 5'-triphosphate", "adenosine-5'-triphosphate", 'adenosine-triphosphate', 'adenylpyrophosphate', 'atp']

    self.assertEqual(metabolite.annotations, expected_annotations)
    self.assertEqual(metabolite.names, expected_names)
    self.assertEqual(metabolite._inchi_string, "InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-3/t4-,6-,7-,10-/m1/s1")


class Test_annotateFull(unittest.TestCase):
  def testBiggAnnotate(self):
    metabolite: MeMoMetabolite = MeMoMetabolite()
    metabolite.set_id("13dpg")
    metabolite.annotations = {'bigg.metabolite': ['13dpg']}
    ret = handleMetabolites([metabolite], db_name = "BiGG", allow_missing_dbs = False)
    expected_annotations = {'bigg.metabolite': ['13dpg'], 'reactome': ['R-ALL-29800'], 'kegg.compound': ['C00236'], 'chebi': ['CHEBI:11881', 'CHEBI:16001', 'CHEBI:1658', 'CHEBI:20189', 'CHEBI:57604'], 'hmdb': ['HMDB62758'], 'inchikey': ['LJQLQCAXBUHEAZ-UWTATZPHSA-J'], 'biocyc': ['META:DPG'], 'metanetx.chemical': ['MNXM261'], 'seed.compound': ['cpd00203']}
    self.assertEqual(metabolite.annotations, expected_annotations)
    self.assertEqual(metabolite.names, ["3-Phospho-D-glyceroyl phosphate"])
    self.assertEqual(ret, AnnotationResult(1, 1, 1))


  def testVMHAnnotate(self):
    metabolite: MeMoMetabolite = MeMoMetabolite()
    metabolite.set_id("10fthf")
    metabolite.set_annotations({'vmhmetabolite': ['10fthf']}, source = "VMH")
    ret = handleMetabolites([metabolite], db_name = "VMH", allow_missing_dbs = False)
    expected_annotations =  {'vmhmetabolite': ['10fthf'], 'bigg.metabolite': ['10fthf'], 'kegg.compound': ['C00234'], 'chemspider': ['109092'], 'chebi': ['15637'], 'biocyc': ['10-FORMYL-THF'], 'foodb.compound': ['FDB022345'], 'hmdb': ['HMDB0000972'], 'metanetx': ['MNXM237'], 'metlin': ['5912'], 'pubchem.compound': ['122347'], 'seed.compound': ['cpd00201'], 'cas': ['2800-34-2']}
    self.assertEqual(metabolite.annotations, expected_annotations)
    print(metabolite.names)
    self.assertEqual(metabolite.names, ['(2S)-2-[(4-{N-[(2-amino-4-oxo-1,4,5,6,7,8-hexahydropteridin-6-yl)methyl]formamido}phenyl)formamido]pentanedioic acid', '10-Formyltetrahydrofolate'])
    self.assertEqual(ret, AnnotationResult(1, 1, 1))


  def testSEEDAnnotate(self):
    metabolite: MeMoMetabolite = MeMoMetabolite()
    metabolite.set_id('cpd00052')
    metabolite.set_annotations({'seed.compound': ['cpd00052']}, source = "test")
    ret = handleMetabolites([metabolite], db_name = "ModelSeed", allow_missing_dbs = False)
    self.assertEqual(ret, AnnotationResult(1, 1, 1))

    self.assertEqual(metabolite.names, sorted(['ctp', "cytidine 5'-triphosphate", 'cytidine triphosphate', "cytidine-5'-triphosphate", 'cytidine-triphosphate']))
    self.assertEqual(metabolite.annotations, {'seed.compound': ['cpd00052'], 'bigg.metabolite': ['ctp'], 'kegg.compound': ['C00063'], 'metacyc.compound': ['CTP']})
    self.assertEqual(metabolite._inchi_string, "InChI=1S/C9H16N3O14P3/c10-5-1-2-12(9(15)11-5)8-7(14)6(13)4(24-8)3-23-28(19,20)26-29(21,22)25-27(16,17)18/h1-2,4,6-8,13-14H,3H2,(H,19,20)(H,21,22)(H2,10,11,15)(H2,16,17,18)/p-3/t4-,6-,7-,8-/m1/s1")


  def testChEBIAnnotate(self):
    metabolite: MeMoMetabolite = MeMoMetabolite()
    metabolite.set_id('117228')
    metabolite.set_annotations({'chebi': ['117228']}, source = "test")
    ret = handleMetabolites([metabolite], db_name ="ChEBI", allow_missing_dbs = False)
    expected_inchi = "InChI=1S/C20H27N3O2/c1-4-20(3)13-14-9-6-7-10-15(14)17-16(20)18(25)23(5-2)19(22-17)21-11-8-12-24/h6-7,9-10,24H,4-5,8,11-13H2,1-3H3,(H,21,22)"
    self.assertEqual(ret, AnnotationResult(1, 0, 0))
    self.assertEqual(metabolite._inchi_string, expected_inchi)

  def testHMDBAnnotate(self):
    metabolite: MeMoMetabolite = MeMoMetabolite()
    metabolite.set_id("10fthf")
    metabolite.set_annotations({'HMDB': ['HMDB0000972']}, source = "test")
    expected_annotations = {'HMDB': ['HMDB0000972'], 'hmdb': ['HMDB0000972'], 'chemspider': ['109092'], 'metlin': ['5912'], 'food.compound': ['FDB030256'], 'pubchem.compound': ['122347'], 'chebi': ['15637'], 'kegg.compound': ['C00234'], 'bigg.metabolite': ['34337'], 'vmhmetabolite': ['10FTHF']}
    ret = handleMetabolites([metabolite], db_name = "HMDB", allow_missing_dbs = False)
    self.assertEqual(metabolite.annotations, expected_annotations)
    print(metabolite.names)
    self.assertEqual(metabolite.names, ['(2S)-2-[(4-{N-[(4-hydroxy-2-imino-5,6,7,8-tetrahydro-1H-pteridin-6-yl)methyl]formamido}phenyl)formamido]pentanedioic acid', '10-Formyltetrahydrofolate'])
    self.assertEqual(ret, AnnotationResult(1, 1, 1))


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
    self.assertEqual(extracted1, expected1)

    expected3 = ({'AraCyc': ['ATP'], 'BiGG': ['atp'], 'BrachyCyc': ['ATP'], 'KEGG': ['C00002'], 'MetaCyc': ['ATP']}, ['ATP', "Adenosine 5'-triphosphate", "adenosine-5'-triphosphate", 'adenosine-triphosphate', 'adenylpyrophosphate'])
    extracted3 = extractModelSEEDAnnotationsFromAlias(aliases3)
    self.assertEqual(extracted3, expected3)
