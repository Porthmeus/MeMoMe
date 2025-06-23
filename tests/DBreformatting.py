# Porthmeus
# 22.05.25 

import unittest
from src.dbhandling.reformatAux import *
import src.dbhandling.reformatVMH as VMH
import src.dbhandling.reformatHMDB as HMDB
from pathlib import Path
import pandas as pd
import os
import shutil
from io import StringIO
import sys

class Test_DBreformatting(unittest.TestCase):
    # The directory of this file
    #this_directory = Path("tests")
    this_directory = Path(__file__).parent
    dat = this_directory.joinpath("dat")

    def setUp(self):
        self.test_dir = "TestCaseDir"
        os.makedirs(self.test_dir, exist_ok = True)

    def tearDown(self):
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

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

        sys.stdout = sys.__stdout__  # Reset redirect

        self.assertIn("Warning: Column 'Z' not found", captured_output.getvalue())
        self.assertListEqual(list(result.columns), ['X', 'B', 'C'])

    def test_no_rename(self):
        result = HMDB.rename_columns_safe(self.df.copy(), {'X': 'Y'})
        self.assertListEqual(list(result.columns), ['A', 'B', 'C'])


if __name__ == '__main__':
    unittest.main()
