import unittest
from contextlib import redirect_stdout
from pathlib import Path
from copy import deepcopy

import cobra
import pandas as pd

from src.ModelMerger2 import ModelMerger
from src.MeMoModel import MeMoModel
from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix

class Test_ModelMerger(unittest.TestCase):
    
    def test_translate(self):
        cmod1 = cobra.io.load_model("textbook")
        mod1 = MeMoModel().fromModel(cmod1)

        cmod2 = cobra.io.load_model("textbook")
        for met in cmod2.metabolites:
            met.id = "M2_" + met.id.replace("_e","_e0")
            if met.compartment == "e":
                met.compartment = "e0"
        for rxn in cmod2.exchanges:
            rxn.id =  rxn.id.replace("EX_","EX_M2_").replace("_e","_e0")
        cmod2._compartments.update({"e0":"external2"})

        mod2 = MeMoModel().fromModel(cmod2)
        mod1.annotated = True
        mod2.annotated = True
        match_table = pd.DataFrame({"met_id1":sorted(list(set([handle_metabolites_prefix_suffix(x.id) for x in cmod1.metabolites]))),
                                    "met_id2":sorted(list(set([handle_metabolites_prefix_suffix(x.id) for x in cmod2.metabolites])))})
        modMerger = ModelMerger(mod1, mod2, match_table)
        modMerger.translate()
        modMerger.merge_models_simple()
#        mergedMod = modMerger.merged_model.copy()
#        mergedMod.notes
#        [x.id for x in mergedMod.reactions if x.objective_coefficient == 1]
#        mergedMod.optimize()
#        for x in mergedMod.reactions:
#            if x.lower_bound > 0:
#                x.lower_bound = 0
#
#        [x.id for x in mergedMod.reactions if x.lower_bound >0]
        self.assertEqual(modMerger.model1.compartments, {"c":"cytosol", "e":"extracellular", "t":"translation"})
        self.assertTrue(all([rxn.compartments == {"e"} for rxn in modMerger.model1.exchanges]))
        self.assertTrue(all([rxn.compartments == {"t","e"} for rxn in modMerger.model1.reactions if rxn.id.startswith("TR_")]))

        self.assertEqual(modMerger.model2.compartments, {"c":"cytosol", "e":"extracellular", "t":"translation"})
        self.assertTrue(all([rxn.compartments == {"e"} for rxn in modMerger.model2.exchanges]))
        self.assertTrue(all([rxn.compartments == {"t","e"} for rxn in modMerger.model2.reactions if rxn.id.startswith("TR_")]))
        self.assertTrue(round(cmod1.optimize().objective_value,3) == round(tmod1.optimize().objective_value,3))
        self.assertTrue(round(cmod2.optimize().objective_value,3) == round(tmod2.optimize().objective_value,3))



#def is_tr_reaction(reaction: cobra.Reaction) -> bool:
#    translation_mets = [
#        met
#        for met in reaction.metabolites
#        if met.compartment == ModelMerger.TRANSLATION_COMPARTMENT
#    ]
#    if len(translation_mets) != 1:
#        return False
#    non_translation = [
#        met
#        for met in reaction.metabolites
#        if met.compartment != ModelMerger.TRANSLATION_COMPARTMENT
#    ]
#    return len(non_translation) == 1
#
#
#def tr_reaction_prefix(reaction_id: str) -> str | None:
#    if reaction_id.startswith("TR_model1_"):
#        return "model1"
#    if reaction_id.startswith("TR_model2_"):
#        return "model2"
#    return None
#
#
#class Test_ModelMerger(unittest.TestCase):
#    this_directory = Path(__file__).parent
#    dat = this_directory / "dat"
#    target_ori = MeMoModel.fromPath(dat / "tiny_ecoli_keep_inchi.xml")
#    source_ori = MeMoModel.fromPath(dat / "tiny_myb11.xml")
#    #target_ori.annotated = True
#    #source_ori.annotated = True
#    target_ori.annotate(allow_missing_dbs=True)
#    source_ori.annotate(allow_missing_dbs=True)
#
#    def setUp(self):
#        self.target = self.target_ori.copy()
#        self.source = self.source_ori.copy()
#
#    def test_translate_namespace_creates_translation_compartment(self):
#       # self.target.annotate(allow_missing_dbs=True)
#       # self.source.annotate(allow_missing_dbs=True)
#        matches = self.target.match(self.source, keepAllMatches=True)
#        merger = ModelMerger(self.target, self.source, matches)
#        merged_model = merger.translate_namespace()
#        self.assertIn(ModelMerger.TRANSLATION_COMPARTMENT, merged_model.compartments)
#
#        translation_reactions = [
#            rxn for rxn in merged_model.reactions if is_tr_reaction(rxn)
#        ]
#        # TODO: Check that there are 2 translation reactions for each matched metabolite plus one translation reaction for each unmatched metabolite in source and target
#
#
#        translation_met_ids = {met.id for met in merged_model.metabolites if met.compartment == ModelMerger.TRANSLATION_COMPARTMENT}
#        for expected_met in ("h2o_t", "o2_t", "pi_t", "co2_t"):
#            self.assertIn(expected_met, translation_met_ids)
#            self.assertIn(f"EX_{expected_met}", [rxn.id for rxn in merged_model.reactions])
#
#        unnamed = [met for met in merged_model.metabolites if met.compartment == ModelMerger.TRANSLATION_COMPARTMENT and met.name is None]
#        if unnamed:
#            print("Translation metabolites with missing names:")
#            for met in unnamed:
#                base_id = met.id[:-2] if met.id.endswith("_t") else met.id
#                model1_candidates = [
#                    m.id
#                    for m in merged_model.metabolites
#                    if m.id.startswith("model1_")
#                    and handle_metabolites_prefix_suffix(m.id[len("model1_") :]) == base_id
#                ]
#                model2_candidates = [
#                    m.id
#                    for m in merged_model.metabolites
#                    if m.id.startswith("model2_")
#                    and handle_metabolites_prefix_suffix(m.id[len("model2_") :]) == base_id
#                ]
#                print(
#                    f"{met.id}: model1={model1_candidates} model2={model2_candidates}"
#                )
#
#
#    def test_translate_namespace_with_match_output(self):
#
#        #self.target.annotate(allow_missing_dbs=True)
#        #self.source.annotate(allow_missing_dbs=True)
#        matches = self.target.match(self.source, keepAllMatches=True)
#
#        expected_pairs = [
#            ("h2o", "cpd00001"),
#            ("o2", "cpd00007"),
#            ("pi", "cpd00009"),
#            ("co2", "cpd00011"),
#        ]
#
#        for target_id, source_id in expected_pairs:
#            pair = matches[
#                (matches["met_id1"] == target_id)
#                & (matches["met_id2"] == source_id)
#            ]
#            self.assertFalse(pair.empty, f"Expected mapping {target_id}<-{source_id} not found")
#        ## assert that no other matches exist
#        #self.assertEqual(len(matches), len(expected_pairs), "Unexpected extra matches found")
#        merger = ModelMerger(self.target, self.source, matches)
#        target_exchange_map = merger._build_exchange_map("model1_")
#        source_exchange_map = merger._build_exchange_map("model2_")
#        expected_matches = merger.matches[
#            merger.matches["target_namespace"].isin(target_exchange_map)
#            & merger.matches["source_namespace"].isin(source_exchange_map)
#        ]
#
#        # select matche which pass the matching score thresholds
#        sel = (
#                (expected_matches["inchi_score"] == 1.0) | 
#                ((expected_matches["Name_score"] >= 0.9) & 
#                 (expected_matches["DB_score"] >= 0.5))
#                )
#        expected_matches = expected_matches.loc[sel]
#       # if "Name_score" in expected_matches.columns:
#       #     expected_matches = expected_matches.loc[expected_matches["Name_score"] >= 0.9]
#       # else:
#       #     self.fail("Name_score not present in matches table columns")
#        expected_matches = (
#            expected_matches.sort_values(
#                by=["total_score", "target_namespace", "source_namespace"],
#                ascending=[False, True, True],
#            )
#            .drop_duplicates(subset=["target_namespace"], keep="first")
#        )
#        expected_targets = set(expected_matches["target_namespace"])
#        merged_model = merger.translate_namespace()
#
#        #print metabolites in translation compartment
#        translation_mets = [met.id for met in merged_model.metabolites if met.compartment == ModelMerger.TRANSLATION_COMPARTMENT]
#        print("Metabolites in translation compartment:")
#        print(translation_mets)
#
#        #print TR_ reactions
#        translation_reactions = [
#            rxn.id for rxn in merged_model.reactions if is_tr_reaction(rxn)
#        ]
#        print("Translation reactions:")
#        print(translation_reactions)
#
#        # print EX_ reactions in translation compartment and their bounds
#        ex_translation_reactions = [
#            rxn for rxn in merged_model.reactions
#            if rxn.id.startswith("EX_") and any(met.compartment == ModelMerger.TRANSLATION_COMPARTMENT for met in rxn.metabolites)
#        ]
#        missing_names = [rxn.id for rxn in ex_translation_reactions if not rxn.name]
#        self.assertFalse(
#            missing_names,
#            f"Missing names for translation exchange reactions: {missing_names}",
#        )
#        print("Exchange reactions in translation compartment:")
#        for rxn in ex_translation_reactions:
#            print(f"{rxn.id}: bounds ({rxn.lower_bound}, {rxn.upper_bound})")       
#            print(f"  stoichiometry: {rxn.metabolites}")
#
#        tr_reactions = [
#            rxn for rxn in merged_model.reactions if is_tr_reaction(rxn)
#        ]
#        print("Translation reactions stoichiometry:")
#        for rxn in tr_reactions:
#            print(f"{rxn.id}: bounds ({rxn.lower_bound}, {rxn.upper_bound})")
#            print(f"  stoichiometry: {rxn.metabolites}")
#
#        translation_mets = [
#            met for met in merged_model.metabolites if met.compartment == ModelMerger.TRANSLATION_COMPARTMENT
#        ]
#        connected_targets = set()
#        for met in translation_mets:
#            base_id = handle_metabolites_prefix_suffix(met.id)
#            source_met = merger._find_source_metabolite(base_id)
#            self.assertIsNotNone(source_met, f"Missing source metabolite for {met.id}")
#            if source_met is not None:
#                if source_met.name is not None:
#                    self.assertEqual(met.name, source_met.name)
#                if source_met.formula is not None:
#                    self.assertEqual(met.formula, source_met.formula)
#                if source_met.charge is not None:
#                    self.assertEqual(met.charge, source_met.charge)
#                if source_met.annotation:
#                    self.assertEqual(met.annotation, source_met.annotation)
#
#            prefixes = set()
#            for rxn in met.reactions:
#                if not is_tr_reaction(rxn):
#                    continue
#                prefix = tr_reaction_prefix(rxn.id)
#                if prefix:
#                    prefixes.add(prefix)
#            if prefixes == {"model1", "model2"}:
#                connected_targets.add(base_id)
#        self.assertSetEqual(
#            connected_targets,
#            expected_targets,
#            "TR_ reactions connect both models only for high-confidence matches",
#        )
#
#        output_dir = self.dat / "output"
#        output_dir.mkdir(parents=True, exist_ok=True)
#        output_path = output_dir / "automatically_merged_tiny.xml"
#        cobra.io.write_sbml_model(merged_model, output_path)
#
#
#if __name__ == "__main__":
#    unittest.main()
#
#
#class Test_ModelMerger_Slow(unittest.TestCase):
#    def test_merge_recon3d_gapseq_models(self):
#        this_directory = Path(__file__).parent
#        dat = this_directory / "dat" / "manually_merged_models" / "gapseq_recon3D"
#        target_path = dat / "M1_recon3D_301_modified.xml"
#        source_path = dat / "M2_bacterial_model.xml"
#        output_dir = dat / "output"
#        output_dir.mkdir(parents=True, exist_ok=True)
#        log_path = output_dir / "automatically_merged_metamodel_prints.txt"
#        output_path = output_dir / "automatically_merged_metamodel.xml"
#
#        with log_path.open("w") as log_file, redirect_stdout(log_file):
#            target_model = cobra.io.read_sbml_model(target_path)
#            source_model = cobra.io.read_sbml_model(source_path)
#            target = MeMoModel.fromModel(target_model)
#            source = MeMoModel.fromModel(source_model)
#            target.annotate(allow_missing_dbs=True)
#            source.annotate(allow_missing_dbs=True)
#            matches = target.match(source, keepAllMatches=True)
#
#            merger = ModelMerger(target, source, matches)
#            merged_model = merger.translate_namespace()
#
#            cobra.io.write_sbml_model(merged_model, output_path)
#
#            # print metabolites in translation compartment
#            translation_mets = [
#                met.id
#                for met in merged_model.metabolites
#                if met.compartment == ModelMerger.TRANSLATION_COMPARTMENT
#            ]
#            print("Metabolites in translation compartment:")
#            print(translation_mets)
#
#            # print TR_ reactions
#            translation_reactions = [
#                rxn.id
#                for rxn in merged_model.reactions
#                if is_tr_reaction(rxn)
#            ]
#            print("Translation reactions:")
#            print(translation_reactions)
#
#            # print EX_ reactions in translation compartment and their bounds
#            ex_translation_reactions = [
#                rxn
#                for rxn in merged_model.reactions
#                if rxn.id.startswith("EX_")
#                and any(
#                    met.compartment == ModelMerger.TRANSLATION_COMPARTMENT
#                    for met in rxn.metabolites
#                )
#            ]
#            print("Exchange reactions in translation compartment:")
#            for rxn in ex_translation_reactions:
#                print(f"{rxn.id}: bounds ({rxn.lower_bound}, {rxn.upper_bound})")
#                print(f"  stoichiometry: {rxn.metabolites}")
#
#            tr_reactions = [
#                rxn for rxn in merged_model.reactions if is_tr_reaction(rxn)
#            ]
#            print("Translation reactions stoichiometry:")
#            for rxn in tr_reactions:
#                print(f"{rxn.id}: bounds ({rxn.lower_bound}, {rxn.upper_bound})")
#                print(f"  stoichiometry: {rxn.metabolites}")
#
#        self.assertGreater(len(merged_model.reactions), 0)
