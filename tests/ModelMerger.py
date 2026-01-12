import unittest
from pathlib import Path

import cobra
import pandas as pd

from src.ModelMerger import ModelMerger
from src.MeMoModel import MeMoModel
from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix


class Test_ModelMerger(unittest.TestCase):
    this_directory = Path(__file__).parent
    dat = this_directory / "dat"

    def setUp(self):
        self.target = cobra.io.read_sbml_model(str(self.dat / "tiny_ecoli_keep_inchi.xml"))
        self.source = cobra.io.read_sbml_model(str(self.dat / "tiny_myb11.xml"))
        self.matches = pd.DataFrame(
            {
                "met_id1": ["h2o", "o2", "pi", "co2"],
                "met_id2": ["cpd00001", "cpd00007", "cpd00009", "cpd00011"],
                "matching_score": [1.0, 1.0, 1.0, 1.0],
                "total_score": [1.0, 1.0, 1.0, 1.0],
            }
        )

    def test_translate_namespace_creates_translation_compartment(self):
        merger = ModelMerger(self.target, self.source, self.matches)
        merged_model = merger.translate_namespace()
        self.assertIn("t", merged_model.compartments)

        translation_reactions = [rxn for rxn in merged_model.reactions if rxn.id.startswith("TR_")]
        self.assertGreater(len(translation_reactions), 0)

        for met in merged_model.metabolites:
            if met.compartment == "t":
                self.assertFalse(met.id.startswith(("model1_", "model2_")))

        translation_met_ids = {met.id for met in merged_model.metabolites if met.compartment == "t"}
        for expected_met in ("h2o_t", "o2_t", "pi_t", "co2_t"):
            self.assertIn(expected_met, translation_met_ids)
            self.assertIn(f"EX_{expected_met}", [rxn.id for rxn in merged_model.reactions])

    def test_translate_namespace_with_match_output(self):
        target_model = MeMoModel.fromPath(self.dat / "tiny_ecoli_keep_inchi.xml")
        source_model = MeMoModel.fromPath(self.dat / "tiny_myb11.xml")
        target_model.annotate(allow_missing_dbs=True)
        source_model.annotate(allow_missing_dbs=True)
        matches = target_model.match(source_model, keepAllMatches=True)

        expected_pairs = [
            ("h2o", "cpd00001"),
            ("o2", "cpd00007"),
            ("pi", "cpd00009"),
            ("co2", "cpd00011"),
        ]
        matches_check = matches.copy()
        matches_check["target_canon"] = matches_check["met_id1"].apply(handle_metabolites_prefix_suffix)
        matches_check["source_canon"] = matches_check["met_id2"].apply(handle_metabolites_prefix_suffix)
        for target_id, source_id in expected_pairs:
            pair = matches_check[
                (matches_check["target_canon"] == target_id)
                & (matches_check["source_canon"] == source_id)
            ]
            self.assertFalse(pair.empty, f"Expected mapping {target_id}<-{source_id} not found")

        merger = ModelMerger(target_model.cobra_model, source_model.cobra_model, matches)
        merged_model = merger.translate_namespace()

        translation_met_ids = {met.id for met in merged_model.metabolites if met.compartment == "t"}
        for expected_met in ("h2o_t", "o2_t", "pi_t", "co2_t"):
            self.assertIn(expected_met, translation_met_ids)
            self.assertIn(f"EX_{expected_met}", [rxn.id for rxn in merged_model.reactions])

        for target_id, source_id in expected_pairs:
            translation_id = f"{target_id}_t"
            translation_met = merged_model.metabolites.get_by_id(translation_id)
            tr_id = f"TR_{translation_id}"
            self.assertIn(tr_id, merged_model.reactions)
            tr_rxn = merged_model.reactions.get_by_id(tr_id)
            self.assertTrue(tr_rxn.id.startswith("TR_"))
            self.assertIn(translation_met, tr_rxn.metabolites)

        target_exchange_bounds = {}
        for rxn in target_model.cobra_model.exchanges:
            if len(rxn.metabolites) != 1:
                continue
            met = next(iter(rxn.metabolites))
            base_id = handle_metabolites_prefix_suffix(met.id)
            target_exchange_bounds.setdefault(base_id, (rxn.lower_bound, rxn.upper_bound))

        bounds_checked = None
        for met in merged_model.metabolites:
            if met.compartment != "t" or not met.id.endswith("_t"):
                continue
            tr_id = f"TR_{met.id}"
            ex_id = f"EX_{met.id}"
            if tr_id not in merged_model.reactions or ex_id not in merged_model.reactions:
                continue
            base_id = met.id[:-2]
            if base_id not in target_exchange_bounds:
                continue
            expected_lb, expected_ub = target_exchange_bounds[base_id]
            if expected_lb == -1000 and expected_ub == 1000:
                continue
            ex_rxn = merged_model.reactions.get_by_id(ex_id)
            tr_rxn = merged_model.reactions.get_by_id(tr_id)
            self.assertAlmostEqual(ex_rxn.lower_bound, expected_lb)
            self.assertAlmostEqual(ex_rxn.upper_bound, expected_ub)
            self.assertEqual(tr_rxn.lower_bound, -1000)
            self.assertEqual(tr_rxn.upper_bound, 1000)
            bounds_checked = (base_id, expected_lb, expected_ub)
            break

        self.assertIsNotNone(
            bounds_checked,
            "No translation exchange found with copied bounds; check TR/EX naming or matches",
        )
        print(
            f"ModelMerger: verified {len(expected_pairs)} mappings; "
            f"bounds copied for {bounds_checked[0]} ({bounds_checked[1]}, {bounds_checked[2]})."
        )

    def test_link_shared_metabolites_requires_exact_ids(self):
        merger = ModelMerger(self.target, self.source, self.matches)
        merged_model = merger.translate_namespace()
        tr_id = "TR_h2o"
        if tr_id not in merged_model.reactions:
            model1_candidates = [
                met.id
                for met in merged_model.metabolites
                if met.id.startswith("model1_")
                and handle_metabolites_prefix_suffix(met.id[len("model1_") :]) == "h2o"
            ]
            self.fail(
                f"{tr_id} missing. Expected a linker TR for h2o, but "
                f"only compartmented model1 metabolites were found: {model1_candidates}"
            )

if __name__ == "__main__":
    unittest.main()
