import unittest
from pathlib import Path

import cobra
import pandas as pd

from src.ModelMerger import ModelMerger


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

if __name__ == "__main__":
    unittest.main()
