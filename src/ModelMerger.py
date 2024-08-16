import cobra
import pandas as pd
import re
from src.MeMoModel import MeMoModel



class ModelMerger:

    def __init__(self,
                 memo_model: MeMoModel = None,
                 matches: pd.DataFrame = None) -> None:
        self.memo_model = memo_model
        self.matches = matches

    def convert_exchange_to_translation(self):
        """converts exchange reactions to namespace translation reactions and adds their respective compartment"""
        cobra_model = self.memo_model.cobra_model
        for ex in cobra_model.exchanges:
            if re.match(r"^EX_", ex.id):
                # Replace "EX_" with "TR_" at the start of the string
                new_id = re.sub(r"^EX_", "TR_", ex.id)
                # Update the reaction ID
                ex.id = new_id
                if len(ex.metabolites.keys()) == 1:
                    met = list(ex.metabolites.keys())[0]
                    # create the metabolite in the translation compartment and give it a +1 stoichiometry (production)
                    met_t = met.copy()
                    #  TODO: replace regular expression with function call for the remove suffix function
                    #   (need to modify the handle_metabolites_prefix_suffix_function)
                    met_t.id = re.sub(r"[^_]+$", "t", met.id)  # substitutes the compartment suffix (all
                                                                            # that follows the last underscore) with
                                                                            # the new compartment's symbol
                    met_t.compartment = "t"
                    ex.add_metabolites({met_t: 1.0})
                else:
                    raise ValueError(ex.id + "should contain only one metabolite")
            else:
                raise ValueError("The exchange reaction " + ex.id + "should start with the 'EX_' prefix")

    def translate_metabolites(self, to_translate, matches_df, score_type):
        """
        Translates a list of metabolites from namespace 2 to the best matching metabolite in namespace 1.

        This function takes a list of metabolite IDs from namespace 2 and attempts to translate each
        to the best matching metabolite ID in namespace 1 based on a matching score. The function ensures
        that no two metabolites from the to_translate list are translated to the same ID in namespace 1.
        If two or more metabolites have the same best match, only the one with the highest score
        (or the first alphanumerically in case of a tie) is assigned. Unassigned metabolites will retain
        their original IDs.

        Parameters:
            to_translate (list): A list of metabolite IDs from namespace 2 that need to be translated.
            matches_df (DataFrame): A DataFrame containing the potential matches between namespaces,
                                    with columns ['met_id2', 'met_id1', score_type].

        Returns:
            dict: A dictionary where keys are original metabolite IDs from namespace 2, and values are
                  the corresponding translated IDs from namespace 1 or the original ID if no match is found.
        """

        # Sort the to_translate list alphanumerically
        to_translate = sorted(to_translate)

        # Sort the dataframe by met_id2, score (descending), and then met_id1 (alphanumeric)
        sorted_matches = matches_df.sort_values(by=['met_id2', score_type, 'met_id1'], ascending=[True, False, True])

        # Create a dictionary to store the best matches
        best_matches = {}

        # Iterate over the to_translate list
        for metabolite in to_translate:
            # Filter the matches for the current metabolite
            matches_for_met = sorted_matches[sorted_matches['met_id2'] == metabolite]

            # Assign the best available match
            for _, row in matches_for_met.iterrows():
                match_id = row['met_id1']
                if match_id not in best_matches.values():
                    best_matches[metabolite] = match_id
                    break
            else:
                # If no match is found, keep the original ID
                best_matches[metabolite] = metabolite

        # Return the dictionary of translations
        return best_matches

    def extract_single_produced_metabolite(self, rxn):
        """
        Extracts the single produced metabolite from a reaction. The function assumes that the reaction
        should produce exactly one metabolite. If the reaction produces more than one or none,
        it raises a ValueError.

        Parameters:
            rxn (cobra.core.Reaction): The COBRApy reaction object from which to extract the produced metabolite.

        Returns:
            str: The metabolite ID with any suffix (like "_t") removed.

        Raises:
            ValueError: If the reaction produces zero or more than one metabolite.
        """
        # Select the "produced" (translated) metabolites
        transl_mets = [k for k, v in rxn.metabolites.items() if v == 1.0]

        # Check that exactly one metabolite is produced
        if len(transl_mets) != 1:
            raise ValueError(
                f"Being a namespace translation reaction, {rxn.id} should produce one and only one metabolite, "
                f"but it produces {len(transl_mets)}.")

        # Get the produced metabolite
        met = transl_mets[0]

        whole_met_id = met.id

        # Remove any suffix, such as "_t", from the metabolite ID
        met_id = re.sub(r"_t$", "", whole_met_id)

        return met_id

    def translate_reactions_and_metabolites_ids(self, score_thr, score_type="Name_score"):
        # TODO: decide best match for each pair of metabolites from the 2 models
        # TODO: split the function into smaller functions
        if score_type not in self.matches.columns: # making sure that the chosen score is specified in the matches table
            raise ValueError("Specified score type: " + score_type + "doesn't exist")
        cobra_model = self.memo_model.cobra_model
        reliable_matches = self.matches.loc[self.matches[score_type] >= score_thr]
        to_translate = []
        for rxn in cobra_model.reactions:
            if rxn.id.startswith("TR_"):
                to_translate = to_translate + [self.extract_single_produced_metabolite(rxn)]
        translation_dict = self.translate_metabolites(to_translate, reliable_matches, score_type)
        for met in cobra_model.metabolites:
            if met.compartment == "t":
                whole_met_id = met.id
                # Remove any suffix, such as "_t", from the metabolite ID
                met_id = re.sub(r"_t$", "", whole_met_id)
                met_t_id = translation_dict[met_id] + "_t"
                met.id = met_t_id
                for rxn in met.reactions:
                    if rxn.id.startswith("TR_"):
                        rxn.id = "TR_" + met_t_id

    def translate_namespace(self):
        self.create_translation_compartment()
        self.convert_exchange_to_translation()
        self.translate_reactions_and_metabolites_ids()
