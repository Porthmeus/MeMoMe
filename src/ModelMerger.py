import cobra
import pandas as pd
import re
from src.MeMoModel import MeMoModel
from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix
#from cobra.core import DictList

class ModelMerger:
    """
    This class aims to facilitate the merging of a model to another model that is based on a different namespace.
    It takes a model and a table of matches between the input model and the target namespace as inputs.
    It doesn't perform a full merging, but creates a new "namespace translation compartment" (denoted by the "_t" suffix
    on metabolite ids) that contains all the exchange metabolites of the input model that could be translated to the
    target namespace. The "t" compartment is connected to the original exchange compartment through a set of "translation
     reactions" (denoted by the "TR_" prefix on reaction ids). Each "TR_" reaction connects one
     external metabolite to its corresponding "_t" metabolite. The medium constraints are then moved from the external
     to the translation compartments
    """
    def __init__(self,
                 memo_model: MeMoModel = None,
                 matches: pd.DataFrame = None) -> None:
        self.memo_model = memo_model
        if matches is None:
            raise ValueError("A DataFrame of matches is expected. Received None as input")
        self.matches = matches.copy()
        self.matches = self.matches.rename(columns={'met_id2': 'source_namespace'})
        self.matches = self.matches.rename(columns={'met_id1': 'target_namespace'})

    def convert_exchange_rxns_to_translation_rxns(self):
        """
        Converts exchange reactions to namespace translation reactions and adds their respective compartment.
        For each exchange reaction, the prefix "EX_" gets replaced with "TR_". A new "translated version" of its
        exchanged metabolite will be created. This will be set stoichiometrically as the product of its exchange
        metabolite. This means that the reaction "consumes" an "exchange" metabolite to produce a "translated" one.
        This "translated" metabolite will keep its original namespace for now, and will be truly translated to the
        target namespace at a later step by a dedicated function (self.translate_ids).  The only role of the "translated"
        metabolite in this function is to take part of the network structure on which the actual namespace translation
        will be applied, converting glc__D_t to the target namespace.

        Example:
        The reaction
            ID: EX_glc__D_e
            Reaction: -1 glc__D_e <--> (exchange with the environment)
        gets converted into
            ID: TR_glc__D_e
            Reaction: - 1 glc__D_e <--> + 1 glc__D_t

        Raises:
        -------
        ValueError: If the exchange reaction does not start with the 'EX_' prefix or contains more than one metabolite.
        """

        cobra_model = self.memo_model.cobra_model
        for ex in cobra_model.exchanges:
            if re.match(r"^EX_", ex.id):
                # ensures that, as one should expect by definition, the exchange reaction involves only one metabolite
                if len(ex.metabolites.keys()) == 1:
                    met = list(ex.metabolites.keys())[0]
                    met_t = met.copy()
                    # removes prefix and suffix from the met id
                    met_t.id = handle_metabolites_prefix_suffix(met_t.id)
                    # adds the translation compartment's suffix to the met id
                    met_t.id = met_t.id + "_t"
                    # assigns the new metabolite to the translation ("t") compartment. As verified in the
                    # test_automatic_compartment_creation() test function, if the compartment did not exist yet, the
                    # execution of this command results in the creation of the compartment
                    met_t.compartment = "t"
                    ex.add_metabolites({met_t: 1.0}) # gives the "t" metabolite a +1 stoichiometry (production). So the
                                                    # stoichiometry of the reaction will be consumption of one external
                                                    # metabolite to produce one "translated version" of it
                else:
                    raise ValueError(ex.id + "should contain only one metabolite")
                # Replace "EX_" with "TR_" at the start of the string
                new_id = re.sub(r"^EX_", "TR_", ex.id)
                # Update the reaction ID
                ex.id = new_id
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
            to_translate (list): A list of metabolites that need to be translated.
            matches_df (DataFrame): A DataFrame containing the potential matches between namespaces,
                                    with columns ['source_namespace', 'target_namespace', score_type].
            score_type (string): A string indicating which column from matches_df contains the matching
                                 score that will be used to rank the metabolite matches.
        """
        # # Sort the metabolites list alphanumerically by their ID
        # to_translate = sorted(to_translate, key=lambda metabolite: metabolite.id)

        # Remove "_t" suffix from the metabolite IDs
        for met in to_translate:
            met.id = re.sub(r"_t$", "", met.id)

        sorted_matches = matches_df.sort_values(by=[score_type, 'source_namespace', 'target_namespace'],
                                                ascending=[False, True, True])

        best_matches = dict()

        # Iterate over the to_translate list
        for met_id in sorted_matches["source_namespace"]:
            if met_id in [m.id for m in to_translate]:
                #print([m.id for m in to_translate])
                #print(met_id)
                # Filter the matches for the current metabolite
                matches_for_met = sorted_matches[sorted_matches['source_namespace'] == met_id]
                # Assign the best available match
                for _, row in matches_for_met.iterrows():
                    match_id = row['target_namespace']
                    if (met_id not in best_matches.keys()) and (match_id not in best_matches.values()):
                        best_matches[met_id] = match_id
                        #print(met_id + ": " + match_id + ", " + str(row[score_type]))
                #else:
                #    # If no match is found, keep the original ID
                #    best_matches[met_id] = met_id
                #    print(met_id + ": " + met_id)

        # apply the id translation to the metabolite and re-appends the translation compartment suffix
        for met in to_translate:
            try:
                met.id = best_matches[met.id] + "_t"
            except KeyError:
                met.id = met.id + "_t"

    def translate_ids(self, score_thr, score_type="Name_score"):
        """
        Translates metabolite IDs in self.cobra_model based on a score threshold and updates corresponding
        translation reaction IDs.

        Parameters:
            score_thr (float): Minimum score required for ID translation.
            score_type (str): Score type to filter matches by, must be a valid column in the matches table

        Raises:
            ValueError: If the specified score type is not found in the matches table or if a metabolite has
                        more than one or no "TR_" reaction.
        """
        if score_type not in self.matches.columns: # makes sure that the chosen score is specified in the matches table
            raise ValueError("Specified score type: " + score_type + "doesn't exist")
        cobra_model = self.memo_model.cobra_model
        reliable_matches = self.matches.loc[self.matches[score_type] >= score_thr]
        #  selects the metabolites to be translated
        to_translate = list()
        for met in cobra_model.metabolites:
            # we translate only metabolites that are taking part in exchange reactions:
            # the "t" compartment has been built by the convert_exchange_rxns_to_translation_rxns to contain all of the
            # exchanged metabolites
            if met.compartment == "t":
                to_translate = to_translate + [met]
        # applies the translation to the metabolites
        self.translate_metabolites(to_translate, reliable_matches, score_type)
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
        self.convert_exchange_rxns_to_translation_rxns()
        self.translate_ids(score_thr=0.5)
        self.create_exchanges()
        self.set_rxn_bounds()
