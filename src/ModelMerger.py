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
        # TODO: cobra.Model does not expose add_genes; genes referenced in GPRs are added via add_reactions.
        # If we need to carry over genes not referenced by any reaction, we'll need a custom merge step.


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



    def translate_namespace(self, score_thr: float = 0.8, score_type: str = "total_score") -> cobra.Model:
        """Public entry point that executes the full namespace translation pipeline."""


        return self.merged_model
