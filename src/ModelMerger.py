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
                    met_t.id = re.sub(r"_e0$", "_t", met.id)
                    met_t.compartment = "t"
                    ex.add_metabolites({met_t: 1.0})
                else:
                    raise ValueError(ex.id + "should contain only one metabolite")
            else:
                raise ValueError("The exchange reaction " + ex.id + "should start with the 'EX_' prefix")

    def translate_reactions_and_metabolites_ids(self):
        cobra_model = self.memo_model.cobra_model
        reliable_matches = self.matches.loc[self.matches["Name_score"] > 0.8]
        print(reliable_matches["met_id2"])
        for rxn in cobra_model.reactions:
            if rxn.id.startswith("TR_"):
                met = list(rxn.metabolites.keys())[1]
                whole_met_id = met.id
                print(whole_met_id)
                met_id = re.sub(r"_t$", "", whole_met_id)
                print(met_id)
                new_met_id = reliable_matches.loc[reliable_matches["met_id2"]==met_id]["met_id1"].values
                print(new_met_id)
                if len(new_met_id) > 0:
                    new_met_id = new_met_id[0] + "_t"
                    met.id = new_met_id
    def translate_namespace(self):
        self.create_translation_compartment()
        self.convert_exchange_to_translation()
        self.translate_reactions_and_metabolites_ids()