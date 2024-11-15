#SteFlor 15.11.2024

import cobra
import pandas as pd
import re
from src.MeMoModel import MeMoModel
from src.removeDuplicateMetabolites import (
    findCommonReactions,
    mergeReactions,
    exchangeMetabolite,
    mergeAndRemove,
    removeDuplicateMetabolites
)
from src import download_db

class ModelMerger:
    """
    A class to merge two metabolic models by translating metabolites from one namespace to another
    using a translation compartment, and integrating shared metabolites.
    """

    def __init__(self,
                 model1: cobra.Model,
                 model2: cobra.Model,
                 matches: pd.DataFrame) -> None:
        """
        Initialize the ModelMerger class.

        Parameters:
            model1 (cobra.Model): The first COBRA model (target namespace).
            model2 (cobra.Model): The second COBRA model (source namespace) to be merged into model1.
            matches (pd.DataFrame): DataFrame containing metabolite matches between source and target namespaces.
                                    Must contain columns 'met_id1', 'met_id2', and score columns.
        """
        #download_db.download()
        self.model1 = model1.copy()  # Target model
        self.model2 = model2.copy()  # Source model
        self.matches = matches.copy()
        # Rename columns for consistency
        self.matches = self.matches.rename(columns={'met_id2': 'source_namespace'})
        self.matches = self.matches.rename(columns={'met_id1': 'target_namespace'})
        # Ensure that the matches DataFrame contains necessary columns
        required_columns = {'source_namespace', 'target_namespace'}
        if not required_columns.issubset(self.matches.columns):
            raise ValueError(f"The matches DataFrame must contain columns: {required_columns}")
        # Preprocess models to handle duplicates
        self.preprocess_models()
        # Add prefixes to both models
        self.add_prefix_to_model_ids(self.model1, prefix='model1_')
        self.add_prefix_to_model_ids(self.model2, prefix='model2_')
        # No need to update matches DataFrame with prefixes since we'll remove them during translation
        # Initialize merged model
        self.merged_model = cobra.Model('merged_model')
        # Add metabolites and reactions from both models
        self.merged_model.add_metabolites(self.model1.metabolites)
        self.merged_model.add_reactions(self.model1.reactions)
        self.merged_model.add_metabolites(self.model2.metabolites)
        self.merged_model.add_reactions(self.model2.reactions)
        # Initialize MeMoModel
        self.memo_model = MeMoModel(self.merged_model)

    def preprocess_models(self):
        """
        Preprocesses model1 and model2 to handle duplicate metabolites and reactions.
        """
        self.model1 = self.remove_duplicates_in_model(self.model1)
        self.model2 = self.remove_duplicates_in_model(self.model2)

    def remove_duplicates_in_model(self, model: cobra.Model):
        """
        Removes duplicate metabolites and reactions from the given model.

        Parameters:
            model (cobra.Model): The COBRA model to process.

        Returns:
            cobra.Model: The cleaned model with duplicates removed.
        """
        # Convert to MeMoModel
        memo_model = MeMoModel(cobra_model=model)
        # Remove duplicates
        memo_model, changes = removeDuplicateMetabolites(memo_model)
        # Return the cleaned cobra model
        return memo_model.cobra_model

    def add_prefix_to_model_ids(self, model: cobra.Model, prefix: str):
        """
        Adds a prefix to all metabolite and reaction IDs in a model to prevent conflicts.

        Parameters:
            model (cobra.Model): The model whose IDs will be prefixed.
            prefix (str): The prefix to add to the IDs.
        """
        # Create a mapping from old metabolite objects to new metabolite objects
        old_metabolite_mapping = {}
        for met in model.metabolites:
            old_id = met.id
            met.id = prefix + met.id
            old_metabolite_mapping[old_id] = met

        # Adjust reaction IDs and update metabolite references
        special_prefixes = ['EX_', 'DM_', 'SK_', 'sink_']  # Exchange, Demand, and Sink reactions
        for rxn in model.reactions:
            original_id = rxn.id
            # Check if the reaction ID starts with any of the special prefixes
            special_prefix = next((sp for sp in special_prefixes if original_id.startswith(sp)), None)
            if special_prefix:
                # For special reactions, insert the model prefix after the special prefix
                rxn.id = special_prefix + prefix + original_id[len(special_prefix):]
            else:
                rxn.id = prefix + original_id

            # Update metabolites in reactions
            new_metabolite_dict = {}
            for met, stoich in rxn.metabolites.items():
                # Get the new metabolite object using the old ID
                old_met_id = met.id  # Since met.id has been updated, we need to get the old ID
                if old_met_id in old_metabolite_mapping:
                    new_met = old_metabolite_mapping[old_met_id]
                else:
                    new_met = met  # Should not happen, but just in case
                new_metabolite_dict[new_met] = stoich

            # Update the metabolites in the reaction
            # First, remove all old metabolites
            rxn.subtract_metabolites(rxn.metabolites)
            # Then, add the new metabolites
            rxn.add_metabolites(new_metabolite_dict)

        # Adjust gene IDs
        for gene in model.genes:
            gene.id = prefix + gene.id

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
                            met_t = met.copy()
                            met_t.id = met_t_id
                            met_t.compartment = "t"
                            # Remove model prefix from metabolite name if needed
                            if met_t.name.startswith(model_prefix):
                                met_t.name = met_t.name[len(model_prefix):]
                            cobra_model.add_metabolites([met_t])
                        # Add the metabolite to the translation reaction
                        ex_copy.add_metabolites({met_t: 1.0})
                        # Remove the model prefix from metabolite IDs in the reaction
                        # (Adjust metabolites in ex_copy)
                        new_metabolites = {}
                        for m, stoich in ex_copy.metabolites.items():
                            if m.id.startswith(model_prefix):
                                m.id = m.id.replace(model_prefix, '', 1)
                            new_metabolites[m] = stoich
                        ex_copy.metabolites = new_metabolites
                        cobra_model.add_reactions([ex_copy])
                    else:
                        raise ValueError(f"{ex_copy.id} should contain only one metabolite")
                else:
                    raise ValueError(f"The exchange reaction {ex.id} should start with the 'EX_' prefix")

    def translate_metabolites(self, to_translate, matches_df, score_type):
        """
        Translates metabolites from source namespace to target namespace based on matching scores.

        Parameters:
            to_translate (list): List of metabolites to translate.
            matches_df (DataFrame): DataFrame containing matches between namespaces.
            score_type (str): The score column to use for matching.
        """
        # Remove "_t" suffix from the metabolite IDs and remove model prefixes
        for met in to_translate:
            met.id = re.sub(r"_t$", "", met.id)
            if met.id.startswith('model1_'):
                met.id = met.id.replace('model1_', '', 1)
            elif met.id.startswith('model2_'):
                met.id = met.id.replace('model2_', '', 1)

        # Remove prefixes from matches DataFrame
        matches_df = matches_df.copy()
        matches_df['source_namespace'] = matches_df['source_namespace'].str.replace('model2_', '', regex=False)
        matches_df['target_namespace'] = matches_df['target_namespace'].str.replace('model1_', '', regex=False)

        # Sort matches by score and IDs
        sorted_matches = matches_df.sort_values(by=[score_type, 'source_namespace', 'target_namespace'],
                                                ascending=[False, True, True])

        best_matches = dict()
        used_target_ids = set()

        # Build best matches without duplicates in target namespace
        for _, row in sorted_matches.iterrows():
            source_id = row['source_namespace']
            target_id = row['target_namespace']
            if source_id in [m.id for m in to_translate] and target_id not in used_target_ids:
                best_matches[source_id] = target_id
                used_target_ids.add(target_id)

        # Apply translation to metabolites
        for met in to_translate:
            source_id = met.id
            if source_id in best_matches:
                new_id = best_matches[source_id] + "_t"
            else:
                new_id = source_id + "_t"

            # Check for existing metabolite
            if new_id in self.merged_model.metabolites:
                existing_met = self.merged_model.metabolites.get_by_id(new_id)
                # Replace metabolite in reactions
                for rxn in list(met.reactions):
                    stoich = rxn.metabolites.pop(met)
                    rxn.add_metabolites({existing_met: stoich})
                # Remove old metabolite
                self.merged_model.metabolites.remove(met)
            else:
                met.id = new_id
                met.compartment = 't'

    def translate_ids(self, score_thr, score_type="matching_score"):
        """
        Translates metabolite IDs in the model based on a score threshold and updates corresponding
        translation reaction IDs.

        Parameters:
            score_thr (float): Minimum score required for ID translation.
            score_type (str): Score type to filter matches by, must be a valid column in the matches table.
        """
        if score_type not in self.matches.columns:
            raise ValueError(f"Specified score type: {score_type} doesn't exist in matches DataFrame")
        cobra_model = self.merged_model
        reliable_matches = self.matches.loc[self.matches[score_type] >= score_thr]
        # Select metabolites to translate
        to_translate = [met for met in cobra_model.metabolites if met.compartment == "t"]
        # Apply translation
        self.translate_metabolites(to_translate, reliable_matches, score_type)
        # Translate IDs of the respective TR_ reactions
        for met in to_translate:
            tr_rxns_for_current_met = [rxn for rxn in met.reactions if rxn.id.startswith("TR_")]
            if len(tr_rxns_for_current_met) != 1:
                raise ValueError(f"{met.id} should have one and only one TR_ reaction, it has {len(tr_rxns_for_current_met)} instead")
            tr_rxn = tr_rxns_for_current_met[0]
            tr_rxn.id = "TR_" + met.id

    def create_exchanges(self):
        """
        Creates exchange reactions for metabolites in the translation compartment.
        """
        cobra_model = self.merged_model
        for met in cobra_model.metabolites:
            if met.compartment == "t":
                reaction_id = "EX_" + met.id
                if reaction_id not in [rxn.id for rxn in cobra_model.reactions]:
                    rxn = cobra.Reaction(id=reaction_id, name=met.name + " exchange", lower_bound=-1000, upper_bound=1000)
                    rxn.add_metabolites({met: -1})
                    cobra_model.add_reactions([rxn])

    def set_rxn_bounds(self):
        """
        Sets the bounds of the new exchange reactions based on the original translation reactions,
        and resets the bounds of translation reactions to default.
        """
        cobra_model = self.merged_model
        for met in cobra_model.metabolites:
            if met.compartment == "t":
                ex_rxn_id = "EX_" + met.id
                tr_rxn_id = "TR_" + met.id
                ex_rxn = cobra_model.reactions.get_by_id(ex_rxn_id)
                tr_rxn = cobra_model.reactions.get_by_id(tr_rxn_id)
                # Transfer bounds from translation reaction to exchange reaction
                ex_rxn.lower_bound = tr_rxn.lower_bound
                ex_rxn.upper_bound = tr_rxn.upper_bound
                # Reset translation reaction bounds to default
                tr_rxn.lower_bound = -1000
                tr_rxn.upper_bound = 1000

    def add_translation_reaction(self, met_id, common_met):
        """
        Adds a translation reaction from a metabolite to the common metabolite in the translation compartment.

        Parameters:
            met_id (str): ID of the original metabolite.
            common_met (Metabolite): Common metabolite in the translation compartment.
        """
        # Remove model prefix from met_id
        met_id_no_prefix = met_id.replace('model1_', '', 1).replace('model2_', '', 1)
        reaction_id = 'TR_' + met_id_no_prefix
        # Check if the reaction already exists
        if reaction_id in [rxn.id for rxn in self.merged_model.reactions]:
            # Reaction already exists, no need to add it again
            return
        cobra_model = self.merged_model
        try:
            original_met = cobra_model.metabolites.get_by_id(met_id)
        except KeyError:
            # Metabolite does not exist in the model
            return
        tr_rxn = cobra.Reaction(id=reaction_id)
        tr_rxn.name = f"Translation reaction for {met_id_no_prefix}"
        tr_rxn.add_metabolites({original_met: -1, common_met: 1})
        cobra_model.add_reactions([tr_rxn])

    def link_shared_metabolites(self):
        """
        Links shared metabolites from both models via the translation compartment.
        """
        cobra_model = self.merged_model
        # Remove prefixes from matches DataFrame
        matches_df = self.matches.copy()
        matches_df['target_namespace'] = matches_df['target_namespace'].str.replace('model1_', '', regex=False)
        matches_df['source_namespace'] = matches_df['source_namespace'].str.replace('model2_', '', regex=False)
        # Use the matches DataFrame to identify shared metabolites
        for _, row in matches_df.iterrows():
            met_id_model1 = 'model1_' + row['target_namespace']
            met_id_model2 = 'model2_' + row['source_namespace']
            # Use target namespace metabolite ID (without prefix) as the common ID
            common_met_id = row['target_namespace'] + '_t'
            # Check if the common metabolite already exists
            if common_met_id in cobra_model.metabolites:
                common_met = cobra_model.metabolites.get_by_id(common_met_id)
            else:
                # Create the common metabolite in the translation compartment
                common_met = cobra.Metabolite(id=common_met_id, compartment='t')
                cobra_model.add_metabolites([common_met])
            # Add translation reactions from each original metabolite to the common metabolite
            self.add_translation_reaction(met_id_model1, common_met)
            self.add_translation_reaction(met_id_model2, common_met)

    def translate_namespace(self, score_thr=0.8, score_type="matching_score"):
        """
        Executes the entire namespace translation process, including conversion of exchange reactions,
        translation of IDs, linking shared metabolites, creation of exchanges, and setting reaction bounds.

        Parameters:
            score_thr (float): Minimum score required for ID translation.
            score_type (str): Score type to filter matches by.
        """
        self.convert_exchange_rxns_to_translation_rxns()
        self.translate_ids(score_thr=score_thr, score_type=score_type)
        self.link_shared_metabolites()
        self.create_exchanges()
        self.set_rxn_bounds()
