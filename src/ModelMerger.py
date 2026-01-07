from __future__ import annotations

import re
from typing import Iterable

import cobra
import pandas as pd

from src.MeMoModel import MeMoModel
from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix
from src.removeDuplicateMetabolites import removeDuplicateMetabolites

TRANSLATION_COMPARTMENT = "t"


class ModelMerger:
    """
    merge two metabolic models by translating metabolites from one namespace to another
    using a translation compartment, and integrating shared metabolites.
    """

    def __init__(
        self,
        model1: cobra.Model,
        model2: cobra.Model,
        matches: pd.DataFrame,
    ) -> None:
        """
        Parameters
        ----------
        model1 :
            Target namespace model; exchange IDs from this model drive the final namespace.
        model2 :
            Source namespace model that will be mapped onto model1.
        matches :
            DataFrame produced by ``MeMoModel.match``. Must contain at least the
            ``met_id1``/``met_id2`` columns plus the score column used for filtering.
        """

        self.model1 = model1.copy()
        self.model2 = model2.copy()
        self.matches = matches.copy()
        # Rename columns for clarity
        self.matches = self.matches.rename(columns={'met_id2': 'source_namespace'})
        self.matches = self.matches.rename(columns={'met_id1': 'target_namespace'})
        # Ensure that the matches DataFrame contains necessary columns
        required_columns = {'source_namespace', 'target_namespace'}
        if not required_columns.issubset(self.matches.columns):
            raise ValueError(f"The matches DataFrame must contain columns: {required_columns}")
        # Preprocess models to handle duplicates
        self.preprocess_models()
        self.add_prefix_to_model_ids(self.model1, "model1_")
        self.add_prefix_to_model_ids(self.model2, "model2_")
        # No need to update matches DataFrame with prefixes since we'll remove them during translation
        # Initialize merged model
        self.merged_model = cobra.Model("merged_model")
        # Add metabolites and reactions from both models
        self.merged_model.add_metabolites(self.model1.metabolites)
        self.merged_model.add_reactions(self.model1.reactions)
        self.merged_model.add_metabolites(self.model2.metabolites)
        self.merged_model.add_reactions(self.model2.reactions)

    def preprocess_models(self) -> None:
        """Remove duplicate metabolites/reactions from both inputs."""

        self.model1 = self._remove_duplicates(self.model1)
        self.model2 = self._remove_duplicates(self.model2)

    @staticmethod
    def _remove_duplicates(model: cobra.Model) -> cobra.Model:
        memo_model = MeMoModel(cobra_model=model)
        memo_model, _ = removeDuplicateMetabolites(memo_model)
        return memo_model.cobra_model

    def add_prefix_to_model_ids(self, model: cobra.Model, prefix: str) -> None:
        """
        Adds a prefix to all metabolite and reaction IDs in a model to prevent conflicts.

        Parameters:
            model (cobra.Model): The model whose IDs will be prefixed.
            prefix (str): The prefix to add to the IDs.
        """
        old_metabolite_mapping = {}
        for met in model.metabolites:
            old_id = met.id
            met.id = prefix + met.id
            old_metabolite_mapping[old_id] = met

        special_prefixes = ("EX_", "DM_", "SK_", "sink_")
        for reaction in model.reactions:
            original_id = reaction.id
            special_prefix = next((sp for sp in special_prefixes if original_id.startswith(sp)), None)
            if special_prefix:
                reaction.id = special_prefix + prefix + original_id[len(special_prefix) :]
            else:
                reaction.id = prefix + original_id

            new_metabolite_dict = {}
            for met, stoich in reaction.metabolites.items():
                old_met_id = met.id  # Since met.id has been updated, we need to get the old ID
                if old_met_id in old_metabolite_mapping:
                    new_met = old_metabolite_mapping[old_met_id]
                else:
                    new_met = met
                new_metabolite_dict[new_met] = stoich

            reaction.subtract_metabolites(reaction.metabolites)
            reaction.add_metabolites(new_metabolite_dict)

        for gene in model.genes:
            gene.id = f"{prefix}{gene.id}"

    def convert_exchange_rxns_to_translation_rxns(self):
        """
        Converts exchange reactions from both models to translation reactions and adds corresponding metabolites
        to the translation compartment.
        """
        cobra_model = self.merged_model
        # Combine the models into a list
        models = [self.model1, self.model2]
        for model in models:
            for ex in model.exchanges:
                if ex.id.startswith('EX_'):
                    # Extract the model prefix from the reaction ID
                    # The exchange reaction IDs are in the format 'EX_modelX_met_id'
                    ex_met_id_with_prefix = ex.id[3:]  # Skip 'EX_'
                    # Determine which model we're processing
                    if ex_met_id_with_prefix.startswith('model1_'):
                        model_prefix = 'model1_'
                    elif ex_met_id_with_prefix.startswith('model2_'):
                        model_prefix = 'model2_'
                    else:
                        raise ValueError(f"Unexpected prefix in exchange reaction ID: {ex.id}")
                    # Copy the reaction to avoid modifying the original model
                    ex_copy = ex.copy()
                    # Remove the model prefix from the metabolite ID
                    met_id = ex_met_id_with_prefix.replace(model_prefix, '', 1)
                    # New translation reaction ID without model prefix
                    new_id = 'TR_' + met_id
                    ex_copy.id = new_id
                    if len(ex_copy.metabolites) == 1:
                        met = next(iter(ex_copy.metabolites))
                        # Remove the model prefix from the metabolite ID
                        met_id_no_prefix = met.id.replace(model_prefix, '', 1)
                        # Create the metabolite in the translation compartment without model prefix
                        met_t_id = met_id_no_prefix + '_t'
                        # Check if the metabolite already exists
                        if met_t_id in cobra_model.metabolites:
                            met_t = cobra_model.metabolites.get_by_id(met_t_id)
                        else:
                            met_t = cobra.Metabolite(
                                id=met_t_id,
                                formula=met.formula,
                                name=met.name[len(model_prefix):] if met.name.startswith(model_prefix) else met.name,
                                compartment='t'
                            )
                            cobra_model.add_metabolites([met_t])
                        # Adjust the reaction metabolites
                        # Remove the original metabolite from the reaction
                        stoich = ex_copy.metabolites[met]
                        ex_copy.subtract_metabolites({met: stoich})
                        # Add the new metabolite to the reaction with the same stoichiometry
                        ex_copy.add_metabolites({met_t: stoich})
                        # Add the updated reaction to the merged model
                        cobra_model.add_reactions([ex_copy])
                    else:
                        raise ValueError(f"{ex_copy.id} should contain only one metabolite")
                else:
                    raise ValueError(f"The exchange reaction {ex.id} should start with the 'EX_' prefix")


    @staticmethod
    def _canonicalize_identifier(identifier: str) -> str:
        base = re.sub(r"_t$", "", identifier)
        return handle_metabolites_prefix_suffix(base)

    def translate_metabolites(self, translation_candidates: Iterable[cobra.Metabolite], matches_df: pd.DataFrame, score_type: str) -> None:
        """Rename translation compartment metabolites based on matches."""

        matches_df = matches_df.copy()
        matches_df["source_namespace"] = matches_df["source_namespace"].str.replace("model2_", "", regex=False)
        matches_df["target_namespace"] = matches_df["target_namespace"].str.replace("model1_", "", regex=False)
        matches_df["source_namespace"] = matches_df["source_namespace"].apply(self._canonicalize_identifier)
        matches_df["target_namespace"] = matches_df["target_namespace"].apply(self._canonicalize_identifier)
        sorted_matches = matches_df.sort_values(
            by=[score_type, "source_namespace", "target_namespace"],
            ascending=[False, True, True],
        )

        best_matches: dict[str, str] = {}
        used_targets: set[str] = set()
        translation_ids = {self._canonicalize_identifier(met.id) for met in translation_candidates}

        for _, row in sorted_matches.iterrows():
            source = row["source_namespace"]
            target = row["target_namespace"]
            if source in translation_ids and target not in used_targets:
                best_matches[source] = target
                used_targets.add(target)

        for met in translation_candidates:
            source_id = self._canonicalize_identifier(met.id)
            new_id = f"{best_matches.get(source_id, source_id)}_t"
            if new_id in self.merged_model.metabolites:
                replacement = self.merged_model.metabolites.get_by_id(new_id)
                for rxn in list(met.reactions):
                    stoich = rxn.metabolites[met]
                    rxn.subtract_metabolites({met: stoich})
                    rxn.add_metabolites({replacement: stoich})
                self.merged_model.metabolites.remove(met)
            else:
                met.id = new_id
                met.compartment = TRANSLATION_COMPARTMENT

    def translate_ids(self, score_thr: float, score_type: str = "matching_score") -> None:
        """Rename translation metabolites and their TR reactions."""

        if score_type not in self.matches.columns:
            raise ValueError(f"{score_type!r} not present in matches table columns")

        translation_mets = [met for met in self.merged_model.metabolites if met.compartment == TRANSLATION_COMPARTMENT]
        reliable = self.matches.loc[self.matches[score_type] >= score_thr]
        self.translate_metabolites(translation_mets, reliable, score_type)

        for met in translation_mets:
            tr_rxns = [rxn for rxn in met.reactions if rxn.id.startswith("TR_")]
            if len(tr_rxns) != 1:
                continue
            tr_rxns[0].id = f"TR_{met.id}"

    def create_exchanges(self) -> None:
        """Add an exchange reaction for every translation metabolite if missing."""

        for met in self.merged_model.metabolites:
            if met.compartment != TRANSLATION_COMPARTMENT:
                continue
            ex_id = f"EX_{met.id}"
            if ex_id in self.merged_model.reactions:
                continue
            ex_rxn = cobra.Reaction(ex_id)
            ex_rxn.lower_bound = -1000
            ex_rxn.upper_bound = 1000
            ex_rxn.add_metabolites({met: -1})
            self.merged_model.add_reactions([ex_rxn])

    def set_rxn_bounds(self) -> None:
        """Copy bounds from translation reactions to their corresponding exchanges."""

        for met in self.merged_model.metabolites:
            if met.compartment != TRANSLATION_COMPARTMENT:
                continue
            tr_id = f"TR_{met.id}"
            ex_id = f"EX_{met.id}"
            if tr_id not in self.merged_model.reactions or ex_id not in self.merged_model.reactions:
                continue
            tr_rxn = self.merged_model.reactions.get_by_id(tr_id)
            ex_rxn = self.merged_model.reactions.get_by_id(ex_id)
            ex_rxn.lower_bound = tr_rxn.lower_bound
            ex_rxn.upper_bound = tr_rxn.upper_bound
            tr_rxn.lower_bound = -1000
            tr_rxn.upper_bound = 1000

    def add_translation_reaction(self, prefixed_met_id: str, common_met: cobra.Metabolite) -> None:
        if prefixed_met_id not in self.merged_model.metabolites:
            return
        met = self.merged_model.metabolites.get_by_id(prefixed_met_id)
        base_id = prefixed_met_id.replace("model1_", "", 1).replace("model2_", "", 1)
        rxn_id = f"TR_{base_id}"
        if rxn_id in self.merged_model.reactions:
            return
        tr_rxn = cobra.Reaction(rxn_id)
        tr_rxn.add_metabolites({met: -1, common_met: 1})
        self.merged_model.add_reactions([tr_rxn])

    def link_shared_metabolites(self) -> None:
        matches_df = self.matches.copy()
        matches_df["target_raw"] = matches_df["target_namespace"].str.replace("model1_", "", regex=False)
        matches_df["source_raw"] = matches_df["source_namespace"].str.replace("model2_", "", regex=False)
        matches_df["target_canon"] = matches_df["target_raw"].apply(self._canonicalize_identifier)
        matches_df["source_canon"] = matches_df["source_raw"].apply(self._canonicalize_identifier)
        for _, row in matches_df.iterrows():
            base_id = row["target_canon"]
            common_id = f"{base_id}_t"
            if common_id in self.merged_model.metabolites:
                common_met = self.merged_model.metabolites.get_by_id(common_id)
            else:
                common_met = cobra.Metabolite(common_id, compartment=TRANSLATION_COMPARTMENT)
                self.merged_model.add_metabolites([common_met])
            self.add_translation_reaction(f"model1_{row['target_raw']}", common_met)
            self.add_translation_reaction(f"model2_{row['source_raw']}", common_met)

    def translate_namespace(self, score_thr: float = 0.8, score_type: str = "matching_score") -> cobra.Model:
        """Public entry point that executes the full namespace translation pipeline."""

        self.convert_exchange_rxns_to_translation_rxns()
        self.translate_ids(score_thr=score_thr, score_type=score_type)
        self.link_shared_metabolites()
        self.create_exchanges()
        self.set_rxn_bounds()
        return self.merged_model
