# Porthmeus
# 22.05.25 

import unittest
from src.dbhandling.reformatAux import *
import src.dbhandling.reformatVMH as VMH
import src.dbhandling.reformatHMDB as HMDB
import src.dbhandling.reformatBiGG as BiGG
from pathlib import Path
import pandas as pd
import os
import shutil
from io import StringIO
import sys
import subprocess as sp
from io import TextIOWrapper
from pandas.testing import assert_frame_equal

class Test_DBreformatting(unittest.TestCase):
    # The directory of this file
    #this_directory = Path("tests")
    this_directory = Path(__file__).parent
    dat = this_directory.joinpath("dat")


    def setUp(self):
        # create a directory to test the download function
        self.test_dir = "TestCaseDir"
        os.makedirs(self.test_dir, exist_ok = True)
        # start a https server for testing the download function locally
        self.https_server = sp.Popen(
                ["python3","-m", "http.server", "8888"],
        cwd=str(self.dat.absolute()),
        stdout=sp.PIPE,
        stderr=sp.PIPE
        )

    def tearDown(self):
        # remove the directory which was needed for the download testing
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
        # close the server
        self.https_server.terminate()

    def test_getData(self):
        # test that unallowed dbs lead to an error
        with self.assertRaises(KeyError):
            getData("foobar")
        
        # test if all implemented databases return pd.dataframes
        self.assertEqual(type(getData("VMH")), type(pd.DataFrame()))
        self.assertEqual(type(getData("BiGG")), type(pd.DataFrame()))
        self.assertEqual(type(getData("ModelSeed")), type(pd.DataFrame()))

    def test_writeData(self):
        # test if we can create an empty file
        writeData(pd.DataFrame(), db = "TestCase")
        self.assertTrue(os.path.exists(os.path.join(self.test_dir,"TestCaseFile.txt")))
    
    ### VMH
    def test_reformatVMH(self):
        # tests only individual functions on a shortened data frame
        vmh = getData("VMH")
        vmh = vmh.iloc[0:100,] # reduce the table to the first 100 entries
        # add an entry with known values to evaluate the code
        vmh.loc[len(vmh)] = ["test", # met_id
                             "test", # abbreviation
                             "", #createdDate
                             "", #updatedDate
                             "test metabolite", #fullName
                             "", #description
                             "metabolite for test***test|assertion metabolite", #synonyms
                             "t-metabolite", #iupac
                             "", #neutralFormula
                             "", #chargedFormula
                             "", #charge
                             "",#avgmolweight
                             "",#monoisotopicweight
                             "miriam1", #miriam
                             "biggId1",#biggId
                             "",#lmId
                             "",#ehmnId
                             "",#hepatonetId
                             "keggId1",#keggId
                             "pubchemId1",#pubChemId
                             "cheBlId1",#cheBlId
                             "",#chembl
                             "InChI=1S-Metabolite",#inchiString
                             "",#inchiKey
                             "SMILE-METABOLITE",#smile
                             "hmdb1",#hmdb
                             "metanetx1",#metanetx
                             "seed1",#seed
                             "",#pdmapName
                             "",#reconMap
                             "",#reconMap3
                             "",#golgimap
                             "",#lysosomemap
                             "",#mitochondrionmap
                             "",#nucleusmap
                             "",#peroxisomemap
                             "",#reticulummap
                             "",#food_db
                             "chemspider1",#chemspider
                             "biocyc1",#biocyc
                             "",#wikipedia
                             "drugbank1",#drugbank
                             "knapsack1",#knapsack
                             "",#phenolExplorer
                             "metlin1",#metlin
                             "casRegistry",#casRegistry
                             "",#epa_id
                             "",#echa_id
                             "iuphar_id1",#iuphar_id
                             "",#fda_id
                             "",#mesh_id
                             "",#chodb_id
                             "",#isHuman
                             ""]#isMicrobe
        # expected database format
        db100 = {'vmhmetabolite': ['test'], 'bigg.metabolite': ['biggId1'], 'kegg.compound': ['keggId1'], 'chemspider': ['chemspider1'], 'chebi': ['cheBlId1'], 'biocyc': ['biocyc1'], 'hmdb': ['hmdb1'], 'iuphar.ligand': ['iuphar_id1'], 'metanetx': ['metanetx1'], 'metlin': ['metlin1'], 'pubchem.compound': ['pubchemId1'], 'seed.compound': ['seed1'], 'drugbank': ['drugbank1'], 'knapsack': ['knapsack1'], 'cas': ['casRegistry']}
        # reformat the vmh database and check different things
        dbs_ref = VMH.reformatVMH(vmh=vmh)
        self.assertTrue(all(dbs_ref.columns[0:4] == ['id', 'name', 'inchi', 'DBs']))
        self.assertEqual(dbs_ref.id.iloc[100], "test")
        self.assertEqual(dbs_ref.name.iloc[100], 'test metabolite_|_t-metabolite_|_metabolite for test_|_test|assertion metabolite')
        self.assertEqual(dbs_ref.inchi.iloc[100], "InChI=1S-Metabolite")
        self.assertEqual(str(db100), dbs_ref.DBs.iloc[100])
    
    # BiGG
    def test_reformatBiGG(self):
        # tests only individual functions on a shortened data frame
        bigg = getData("BiGG")
        bigg.columns
        bigg = bigg.iloc[0:100,] # reduce the table to the first 100 entries

        # add a line with known values
        bigg.loc[len(bigg)] = ["bigg_id", # met_id
                               "test", # universal_bigg_id 
                               "test metabolite", # name
                               "", # model_list 
                               "Test:  http://identifiers.org/test/test1", # database_links
                               "bigg_id_c; bigg_id[c]"] # old_bigg_ids 
        # add another line with slight differences, to test the aggregation of universal_bigg_ids
        bigg.loc[len(bigg)] = ["bigg_id", # met_id
                               "u_bigg_id", # universal_bigg_id 
                               "test|metabolite", # name
                               "", # model_list - not important, will be discarded
                               "Test:  http://identifiers.org/test/test2; InChI key https://identifiers.org/inchikey/VHRGRCVQAFMJIZ-UHFFFAOYSA-P", # database_links
                               "bigg_id_c; bigg_id(c)"] # old_bigg_ids 
        # reformat and test
        dbs_ref = BiGG.reformatBiGG(dat = bigg, inchiKey2String_url = "http://localhost:8888/TestInchiKey2String.gz")
        
        db100 = {'test': ['test1','test2']}
        self.assertTrue(all(dbs_ref.columns[0:4] == ['id', 'name', 'inchi', 'DBs']))
        self.assertEqual(dbs_ref.id.loc["test"], "test")
        self.assertEqual(dbs_ref.name.loc["test"], 'test metabolite_|_test|metabolite')
        self.assertEqual(dbs_ref.inchi.loc["test"], "InChI=1S/C5H14N2/c6-4-2-1-3-5-7/h1-7H2/p+2")
        self.assertEqual(str(db100), dbs_ref.DBs.loc["test"])


import unittest
import pandas as pd
import warnings

class TestConcatCols(unittest.TestCase):

    def setUp(self):
        self.row = pd.Series({
            'a': 'foo',
            'b': 'bar',
            'c': '',
            'd': None,
            'e': 'baz|qux'
        })

    def test_basic_concat(self):
        result = HMDB.concatCols(self.row, ['a', 'b'])
        self.assertEqual(result, 'foo|bar')

    def test_ignore_empty_and_none(self):
        result = HMDB.concatCols(self.row, ['a', 'c', 'd', 'b'])
        self.assertEqual(result, 'foo|bar')

    def test_empty_result(self):
        result = HMDB.concatCols(self.row, ['c', 'd'])
        self.assertEqual(result, '')

    def test_separator_replacement_warning(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = HMDB.concatCols(self.row, ['a', 'e'])
            self.assertEqual(result, 'foo|baz*/*qux')
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[0].category, UserWarning))
            self.assertIn("Replacing | in", str(w[0].message))

class TestRenameColumnsSafe(unittest.TestCase):

    def setUp(self):
        self.df = pd.DataFrame({
            'A': [1, 2],
            'B': [3, 4],
            'C': [5, 6]
        })

    def test_successful_rename(self):
        result = HMDB.rename_columns_safe(self.df.copy(), {'A': 'X', 'B': 'Y'})
        self.assertListEqual(list(result.columns), ['X', 'Y', 'C'])

    def test_partial_rename_with_warning(self):
        df_copy = self.df.copy()
        captured_output = StringIO()
        sys.stdout = captured_output  # Redirect stdout

        result = HMDB.rename_columns_safe(df_copy, {'A': 'X', 'Z': 'W'})

        self.assertIn("Warning: Column 'Z' not found", captured_output.getvalue())
        self.assertListEqual(list(result.columns), ['X', 'B', 'C'])

    def test_no_rename(self):
        result = HMDB.rename_columns_safe(self.df.copy(), {'X': 'Y'})
        self.assertListEqual(list(result.columns), ['A', 'B', 'C'])


class TestPrepare(unittest.TestCase):

    def test_prepare_normal_case(self):
        df = pd.DataFrame({
            "name": ["Water", "Ethanol"],
            "synonyms": ["Dihydrogen monoxide", "Ethyl alcohol"],
            "iupac_names": [None, "Ethanol"],
            "traditional_iupac": ["Oxidane", None]
        })
        expected = pd.DataFrame({
            "name": [
                "Water_|_Dihydrogen monoxide_|_Oxidane",
                "Ethanol_|_Ethyl alcohol_|_Ethanol"
            ]
        })
        result = HMDB.prepare(df)
        assert_frame_equal(result, expected)


    def test_prepare_none_input(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = HMDB.prepare(None)
            self.assertIsNone(result)
            self.assertTrue(any("Error during concatenation" in str(warn.message) for warn in w))

class TestGetAnnos(unittest.TestCase):
    this_directory = Path(__file__).parent
    dat_file = this_directory / "dat" / "hmdb_metab_test.xml"
    
    @classmethod
    def setUpClass(cls):
      with open(TestGetAnnos.dat_file) as xml_file:
          cls.df = xml_to_pandas_lazy(xml_file, "metabolite")
          cls.df = HMDB.prepare(cls.df)
          cls.df = HMDB.rename_columns_safe(cls.df, HMDB._keys)


    def test_getAnnosPerEntry(self):
      ret = HMDB.getAnnosPerEntry(self.df, "HMDB0000001")
      self.assertEqual(ret, {'hmdb': ['HMDB0000001'], 'chemspider': ['83153'], 'drugbank': ['DB04151'], 'metlin': ['3741'], 'food.compound': ['FDB093588'], 'pubchem.compound': ['92105'], 'chebi': ['50599'], 'kegg.compound': ['C01152']})

    def test_getAnnos(self):
      ret = pd.DataFrame(HMDB.getAnnos(self.df))
      expected = pd.read_csv(self.this_directory / "dat" / "hmdb_test_reformated.csv", header = None)
      expected = pd.DataFrame(expected.iloc[:, 1])
      ret.reset_index(drop=True, inplace=True)

      ret.columns = range(ret.shape[1])
      expected.columns = range(expected.shape[1])

      assert_frame_equal(ret, expected)
    
    def test_doAll(self):
      df = None
      with open(TestGetAnnos.dat_file) as xml_file:
        df = xml_to_pandas_lazy(xml_file, "metabolite")
      ret = pd.DataFrame(HMDB.do_all(df))
      expected = pd.read_csv(self.this_directory / "dat" / "hmdb_metabolites_reformatted_all.csv", index_col = 0)
      pd.set_option('display.max_rows', None)
      pd.set_option('display.max_columns', None)
      print(ret)
      print(expected)
      assert_frame_equal(ret, expected)





if __name__ == '__main__':
    unittest.main()
