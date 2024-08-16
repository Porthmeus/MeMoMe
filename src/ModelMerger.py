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
                    #  TODO: replace regular expression with function call for the remove suffix function (need to modify the handle_metabolites_prefix_suffix_function)
                    met_t.id = re.sub(r"[^_]+$", "t", met.id)  # substitutes the compartment suffix (all
                                                                            # that follows the last underscore) with
                                                                            # the new compartment's symbol
                    met_t.compartment = "t"
                    ex.add_metabolites({met_t: 1.0})
                else:
                    raise ValueError(ex.id + "should contain only one metabolite")
            else:
                raise ValueError("The exchange reaction " + ex.id + "should start with the 'EX_' prefix")

    def translate_reactions_and_metabolites_ids(self, score_thr, score_type="Name_score"):
        if score_type not in self.matches.columns:  # making sure that the chosen score is specified in the matches table
            raise ValueError("Specified score type: " + score_type + "doesn't exist")
        cobra_model = self.memo_model.cobra_model
        reliable_matches = self.matches.loc[self.matches[score_type] > score_thr]
        for rxn in cobra_model.reactions:
            if rxn.id.startswith("TR_"):
                transl_mets = [k for k, v in rxn.metabolites.items()
                               if v == 1.0]  # select the "produced" (translated) metabolites
                if len(transl_mets) != 1:
                    raise ValueError("Being a namespace translation reaction, " + rxn.id +
                                     "should produce one and only one metabolite, but it produces "
                                     + str(len(transl_mets)))
                met = transl_mets[0]  # I select the one (and only) metabolite produced by the TR_ reaction
                whole_met_id = met.id
                met_id = re.sub(r"_t$", "", whole_met_id)
                new_met_id = reliable_matches.loc[reliable_matches["met_id2"] ==
                                                  met_id].sort_values(by=score_type,
                                                                      ascending=False)["met_id1"].values
                if len(new_met_id) > 0:  # select the first (most reliable mathc) metabolite id from the matches table
                    new_met_id = new_met_id[0] + "_t"
                    met.id = new_met_id

    def translate_namespace(self):
        self.create_translation_compartment()
        self.convert_exchange_to_translation()
        self.translate_reactions_and_metabolites_ids()
