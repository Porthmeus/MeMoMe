#SteFlor 15.11.2024

import unittest
import cobra
import pandas as pd
import os
from src.ModelMerger import ModelMerger
from src.MeMoModel import MeMoModel
from pathlib import Path


class TestModelMerger(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Set up paths to the test models
        # Replace these paths with the actual paths to your test models
        cls.this_directory = Path(__file__).parent
        cls.model1_path = cls.this_directory / 'dat/tiny_ecoli_keep_inchi_withduplicates.xml'
        cls.model2_path = cls.this_directory / 'dat/concoct_1.xml'

        # Load the models
        cls.model1 = cobra.io.read_sbml_model(cls.model1_path)
        cls.model2 = cobra.io.read_sbml_model(cls.model2_path)

        # Load or create the matches DataFrame
        # Replace this with your actual matches DataFrame
        # For demonstration, we'll create a sample DataFrame with placeholder data
        cls.matches = pd.DataFrame({ #TODO: fill with the metabolites that we know are true matches
            'met_id1': [],
            'met_id2': [],
            'matching_score': []
        })

        # Ensure the matches DataFrame has the correct columns
        cls.matches = cls.matches.rename(columns={'met_id1': 'met_id1', 'met_id2': 'met_id2'})

    def test_model_merger_initialization(self):
        """
        Test the initialization of the ModelMerger class.
        """
        # Initialize the ModelMerger
        merger = ModelMerger(self.model1, self.model2, self.matches)

        # Check that the merged model is created
        self.assertIsNotNone(merger.merged_model)

        # Check that prefixes have been added correctly
        for met in merger.model1.metabolites:
            self.assertTrue(met.id.startswith('model1_'))
        for met in merger.model2.metabolites:
            self.assertTrue(met.id.startswith('model2_'))
        for rxn in merger.model1.reactions:
            if rxn.id.startswith(('EX_','DM_', 'SK_', 'sink_')):
                self.assertTrue('model1_' in rxn.id)
            else:
                self.assertTrue(rxn.id.startswith('model1_'))
        for rxn in merger.model2.reactions:
            if rxn.id.startswith(('EX_','DM_', 'SK_', 'sink_')):
                self.assertTrue('model2_' in rxn.id)
            else:
                self.assertTrue(rxn.id.startswith('model2_'))

    #def test_duplicate_removal(self):
    #    """
    #    Test that duplicates are removed during preprocessing.
    #    """
    #    return

    def test_translate_namespace(self):
        """
        Test the translate_namespace method of the ModelMerger class.
        """
        # Initialize the ModelMerger
        merger = ModelMerger(self.model1, self.model2, self.matches)

        # Perform the namespace translation
        merger.translate_namespace(score_thr=0.9, score_type='matching_score')

        # Access the merged model
        merged_model = merger.merged_model

        # Check that the translation compartment exists
        self.assertIn('t', merged_model.compartments)

        # Check that translation reactions have been added
        tr_rxns = [rxn for rxn in merged_model.reactions if rxn.id.startswith('TR_')]
        self.assertGreater(len(tr_rxns), 0)

        # Check that exchange reactions have been created for translation compartment metabolites
        ex_rxns = [rxn for rxn in merged_model.reactions if rxn.id.startswith('EX_')]
        self.assertGreater(len(ex_rxns), 0)
        for rxn in ex_rxns:
            self.assertEqual(rxn.metabolites[next(iter(rxn.metabolites))].compartment, 't')

        # Check that metabolites in the translation compartment do not have model prefixes
        for met in merged_model.metabolites:
            if met.compartment == 't':
                self.assertFalse(met.id.startswith(('model1_', 'model2_')))
        # TODO: implement further checks
        return

    def test_prefixes_in_reactions(self):
        """
        Test that reaction IDs are correctly prefixed, especially for exchange, demand, and sink reactions.
        """
        # Initialize the ModelMerger
        merger = ModelMerger(self.model1, self.model2, self.matches)

        # After merging, check that reactions in the merged model have the correct IDs
        for rxn in merger.merged_model.reactions:
            original_id = rxn.id
            if original_id.startswith('TR_'):
                # Translation reactions should not have model prefixes
                self.assertFalse('model1_' in original_id)
                self.assertFalse('model2_' in original_id)
            elif original_id.startswith('EX_'):
                # Exchange reactions in the translation compartment
                self.assertFalse('model1_' in original_id)
                self.assertFalse('model2_' in original_id)
            else:
                # Other reactions should have model prefixes
                self.assertTrue(original_id.startswith(('model1_', 'model2_')))

    def test_metabolite_ids(self):
        """
        Test that metabolite IDs are correctly handled, especially in the translation compartment.
        """
        # Initialize the ModelMerger
        merger = ModelMerger(self.model1, self.model2, self.matches)
        merger.translate_namespace(score_thr=0.9, score_type='matching_score')

        merged_model = merger.merged_model

        # Check that metabolites in the translation compartment have IDs without model prefixes
        for met in merged_model.metabolites:
            if met.compartment == 't':
                self.assertFalse(met.id.startswith(('model1_', 'model2_')))

        # Check that metabolites in other compartments have the correct prefixes
        for met in merged_model.metabolites:
            if met.compartment != 't':
                self.assertTrue(met.id.startswith(('model1_', 'model2_')))

    def test_gene_ids(self):
        """
        Test that gene IDs are correctly prefixed.
        """
        # Initialize the ModelMerger
        merger = ModelMerger(self.model1, self.model2, self.matches)

        # Check that gene IDs have been prefixed
        for gene in merger.model1.genes:
            self.assertTrue(gene.id.startswith('model1_'))
        for gene in merger.model2.genes:
            self.assertTrue(gene.id.startswith('model2_'))

        # Check genes in the merged model
        for gene in merger.merged_model.genes:
            self.assertTrue(gene.id.startswith(('model1_', 'model2_')))


if __name__ == '__main__':
    unittest.main()
