import cobra
import pandas as pd
import re
from src.MeMoModel import MeMoModel
#from cobra.core import DictList



class ModelMerger:

    def __init__(self,
                 memo_model: MeMoModel = None,
                 matches: pd.DataFrame = None) -> None:
        self.memo_model = memo_model
        self.matches = matches.copy()
        self.matches = self.matches.rename(columns={'met_id2': 'source_namespace'})
        self.matches = self.matches.rename(columns={'met_id1': 'target_namespace'})


    def replace_exchange_rxns_with_translation_rxns(self):
        """copies exchange reactions into new reactions that gets modified to behave as translation reactions. It also
        adds their respective 't' (translation) compartment"""
        cobra_model = self.memo_model.cobra_model
        for ex in cobra_model.exchanges:
            if re.match(r"^EX_", ex.id):
                tr_reac = cobra.Reaction(id=re.sub(r"^EX_", "TR_", ex.id))
                tr_reac.name = ex.name
                tr_reac.lower_bound = ex.lower_bound
                tr_reac.upper_bound = ex.upper_bound
                # copy notes and annotation dictionaries if they exist.
                tr_reac.notes = ex.notes.copy() if ex.notes else {}
                tr_reac.annotation = ex.annotation.copy() if ex.annotation else {}

                # Since this is an exchange reaction, we expect only one metabolite.
                # Iterate over the original reaction’s metabolites (should be one) and clone it.
                if len(ex.metabolites.keys()) == 1:
                    met = list(ex.metabolites.keys())[0]
                    stoich = ex.metabolites[met]
                    # Add the metabolite to the cloned reaction with the same stoichiometry
                    tr_reac.add_metabolites({met: stoich})
                    #remove original exchange reaction
                    cobra_model.reactions.remove(ex)
                    # create the metabolite in the translation compartment and give it the opposite stoichiometry
                    met_t = met.copy()
                    #  TODO: replace regular expression with function call for the remove suffix function
                    #   (need to modify the handle_metabolites_prefix_suffix_function)
                    met_t.id = re.sub(r"[^_]+$", "t", met.id)  # substitutes the compartment suffix (all
                                                                            # that follows the last underscore) with
                                                                            # the new compartment's symbol
                    met_t.compartment = "t"
                    tr_reac.add_metabolites({met_t: -stoich})
                else:
                    raise ValueError(ex.id + "should contain only one metabolite")
            else:
                raise ValueError("The exchange reaction " + ex.id + "should start with the 'EX_' prefix")
            cobra_model.add_reactions([tr_reac])


    def translate_metabolites(self, to_translate, score_type):
        """
        Translates a list of metabolites from namespace 2 to the best matching metabolite in namespace 1.

        This function takes a list of metabolite IDs from namespace 2 and translates each
        to their matching metabolite ID in namespace 1.
        Unassigned (with no match) metabolites will retain their original IDs.

        Parameters:
            to_translate (list): A list of metabolites that need to be translated.
            score_type (string): A string indicating which column from matches_df contains the matching
                                 score that will be used to rank the metabolite matches.
        """

        # Remove "_t" suffix from the metabolite IDs
        for met in to_translate:
            met.id = re.sub(r"_t$", "", met.id)

        # Iterate over the to_translate list to apply the id translation to the metabolite and to re-append the
        # translation compartment suffix
        for met in to_translate:
            try:
                met.id = self.matches[met.id] + "_t"
            except KeyError:
                met.id = met.id + "_t"

    def translate_ids(self, score_type="total_score"):
        """
        Translates metabolite IDs in self.cobra_model based on a score threshold and updates corresponding
        translation reaction IDs.

        Parameters:
            score_type (str): Score type to filter matches by, must be a valid column in the matches table
        Raises:
            ValueError: If the specified score type is not found in the matches table or if a metabolite has
                        more than one or no "TR_" reaction.
        """
        if score_type not in self.matches.columns: # makes sure that the chosen score is specified in the matches table
            raise ValueError("Specified score type: " + score_type + "doesn't exist")
        cobra_model = self.memo_model.cobra_model
        #  selects the metabolites to be translated
        to_translate = list()
        for met in cobra_model.metabolites:
            if met.compartment == "t":
                to_translate = to_translate + [met]
        # applies the translation to the metabolites
        self.translate_metabolites(to_translate, score_type)
        #  translate the ids of the respective TR_ reactions
        for met in cobra_model.metabolites:
            if met.compartment == "t":
                tr_rxns_for_current_met = [rxn for rxn in met.reactions if rxn.id.startswith("TR_")]
                if len(tr_rxns_for_current_met) != 1:
                    raise ValueError(met.id + " should have one and only one TR_ reaction, it has " +
                                     str(len(tr_rxns_for_current_met)) + " instead")
                tr_rxn = tr_rxns_for_current_met[0]
                tr_rxn.id = "TR_" + met.id

    def create_exchanges(self):
        cobra_model = self.memo_model.cobra_model
        for met in cobra_model.metabolites:
            if met.compartment == "t":
                reaction_id = "EX_" + met.id
                rxn = cobra.Reaction(id=reaction_id, name=met.name + " exchange", lower_bound=-1000, upper_bound=1000)
                rxn.add_metabolites({met: -1})
                cobra_model.add_reactions([rxn])

    def set_rxn_bounds(self):
        cobra_model = self.memo_model.cobra_model
        for met in cobra_model.metabolites:
            if met.compartment == "t":
                ex_rxn = cobra_model.reactions.get_by_id("EX_" + met.id)
                tr_rxn = cobra_model.reactions.get_by_id("TR_" + met.id)
                ex_rxn.lower_bound = tr_rxn.lower_bound
                ex_rxn.upper_bound = tr_rxn.upper_bound
                tr_rxn.lower_bound = -1000
                tr_rxn.upper_bound = 1000


    def translate_namespace(self):
        self.replace_exchange_rxns_with_translation_rxns()
        self.translate_ids()
        self.create_exchanges()
        self.set_rxn_bounds()
