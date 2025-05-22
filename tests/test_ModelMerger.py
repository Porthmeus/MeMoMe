import unittest

from cobra.io import write_sbml_model

from src.MeMoModel import MeMoModel
from src.ModelMerger import ModelMerger
from pathlib import Path
import pandas as pd
from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix
from cobra.io import write_sbml_model


class TestModelMerging(unittest.TestCase):
    this_directory = Path(__file__).parent
    dat = this_directory

    def test_translate_namespace(self):
        # Load test models
        model1_path = self.dat.joinpath("dat/tiny_ecoli_keep_inchi.xml")
        model2_path = self.dat.joinpath("manually_merged_models/gapseq_recon3D/M2_bacterial_model.xml")
        model1 = MeMoModel.fromPath(model1_path)
        model2 = MeMoModel.fromPath(model2_path)

        # Hardcoded matches for test simplicity (assumed pre-processed IDs)
        matches = pd.DataFrame({
            'met_id1': ['glc_D'],
            'met_id2': ['cpd00027']
        })

        merger = ModelMerger(model1, model2, matches, fail_on_missing_metabolite=False)
        cobra_model1 = merger.memo_model1.cobra_model

        # Collect original exchange bounds before translation by identifying exchange reactions
        original_bounds = {}
        for source in matches['met_id1']:
            found = False
            for rxn in cobra_model1.reactions:
                if rxn.id.startswith("EX_") and len(rxn.metabolites) == 1:
                    met, stoich = next(iter(rxn.metabolites.items()))
                    # Compare processed ID rather than raw ID suffix
                    if met.compartment == "e" and stoich == -1 and handle_metabolites_prefix_suffix(met.id) == source:
                        original_bounds[source] = (rxn.lower_bound, rxn.upper_bound)
                        found = True
                        break
            if not found:
                self.fail(f"Original exchange reaction for metabolite {source} not found before translation")

        # Perform namespace translation
        merger.translate_namespace()

        # Verify translation reactions and new exchange bounds
        for source, target in zip(matches['met_id1'], matches['met_id2']):
            tr_id = f"TR_{source}"
            # check if translation reaction exists
            try:
                tr_rxn = cobra_model1.reactions.get_by_id(tr_id)
            except KeyError:
                self.fail(f"Missing translation reaction {tr_id}")

            # Find source metabolite without next()
            met_src = None
            for m in cobra_model1.metabolites:
                if handle_metabolites_prefix_suffix(m.id) == source:
                    met_src = m
                    break
            if met_src is None:
                self.fail(f"Source metabolite {source} not found in model1 metabolites")
                        # Find target metabolite without next()
            met_tgt = None
            for m in cobra_model1.metabolites:
                if handle_metabolites_prefix_suffix(m.id) == target:
                    met_tgt = m
                    break
            if met_tgt is None:
                self.fail(f"Target metabolite {target} not found in model1 metabolites")
            stoich = tr_rxn.metabolites
            self.assertEqual(stoich.get(met_src), -1,
                             f"Incorrect stoichiometry for {met_src.id} in {tr_id}")
            self.assertEqual(stoich.get(met_tgt), 1,
                             f"Incorrect stoichiometry for {met_tgt.id} in {tr_id}")


            # Check new exchange reaction bounds
            new_ex_id = f"EX_{target}_e"
            assert new_ex_id in cobra_model1.reactions, f"Missing exchange reaction {new_ex_id}"
            new_ex_rxn = cobra_model1.reactions.get_by_id(new_ex_id)

            lb, ub = original_bounds[source]
            assert new_ex_rxn.lower_bound == lb, \
                f"Lower bound mismatch for {new_ex_id}: expected {lb}, got {new_ex_rxn.lower_bound}"
            assert new_ex_rxn.upper_bound == ub, \
                f"Upper bound mismatch for {new_ex_id}: expected {ub}, got {new_ex_rxn.upper_bound}"

        # Ensure original exchanges removed
        for source in matches['met_id1']:
            for rxn in cobra_model1.reactions:
                if rxn.id.startswith("EX_") and len(rxn.metabolites) == 1:
                    met, stoich = next(iter(rxn.metabolites.items()))
                    if met.compartment == "e" and stoich == -1 and \
                       handle_metabolites_prefix_suffix(met.id) == source:
                        self.fail(f"Original exchange reaction for {source} should have been removed")

        # Non-matching metabolite exchanges remain (EX_ prefix, single e-compartment metabolite with stoich -1)
        non_matching = []
        for rxn in cobra_model1.reactions:
            if rxn.id.startswith("EX_") and len(rxn.metabolites) == 1:
                met, stoich = next(iter(rxn.metabolites.items()))
                if met.compartment == "e" and stoich == -1 and met.id not in matches['met_id1'].tolist():
                    non_matching.append(rxn)
        self.assertTrue(non_matching,
                        "Non-matching metabolite exchange was incorrectly removed")

        print("test_translate_namespace passed")
