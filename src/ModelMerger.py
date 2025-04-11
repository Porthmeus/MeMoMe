import cobra
import pandas as pd
from src.MeMoModel import MeMoModel
from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix

class ModelMerger:
    def __init__(self,
                 memo_model1: MeMoModel = None,
                 memo_model2: MeMoModel = None,
                 matches: pd.DataFrame = None,
                 fail_on_missing_metabolite: bool = True) -> None:
        self.memo_model1 = memo_model1
        self.memo_model2 = memo_model2
        self.matches = matches.copy().rename(columns={'met_id1': 'source_namespace',
                                                      'met_id2': 'target_namespace'})
        self.fail_on_missing_metabolite = fail_on_missing_metabolite

        # Precompute a mapping of processed external metabolite IDs to their cobra.Metabolite in model1.
        self.external_metabolite_map = {}
        cobra_model1 = self.memo_model1.cobra_model
        for met in cobra_model1.metabolites:
            if met.compartment == "e":
                key = handle_metabolites_prefix_suffix(met.id)
                self.external_metabolite_map[key] = met

    def create_translation_reaction(self, met_mod1: cobra.Metabolite, copy_met_mod2: cobra.Metabolite):
        """
        Creates a translation reaction in model1:
          - Copies met_mod2 (from model2) to create copy_met_mod2 and adds it to model1.
          - Creates a reaction with ID "TR_" + <processed met_mod1 id> that converts met_mod1 into the copy of met_mod2.
          - Returns the created reaction and the copy of met_mod2.
        """
        cobra_model1 = self.memo_model1.cobra_model

        # add the copy of met_mod2 to model1.

        cobra_model1.add_metabolites([copy_met_mod2])

        processed_met_mod1_id = handle_metabolites_prefix_suffix(met_mod1.id)
        reaction_id = "TR_" + processed_met_mod1_id
        translation_reaction = cobra.Reaction(reaction_id)
        translation_reaction.add_metabolites({
            met_mod1: -1,      # Consumes the external metabolite in model1.
            copy_met_mod2: 1   # Produces the copied metabolite from model2.
        })
        cobra_model1.add_reactions([translation_reaction])

    def delete_exchange_reaction(self, met_mod1: cobra.Metabolite) -> tuple:
        """
        Deletes any exchange reaction(s) in model1 that exchange met_mod1.
        Assumes an exchange reaction has an ID that starts with "EX_" and involves only met_mod1.
        Returns the lower and upper bounds from the first found exchange reaction,
        or None if none is found.
        """
        cobra_model1 = self.memo_model1.cobra_model

        exchange_rxns = [
            rxn for rxn in cobra_model1.reactions
            if rxn.id.startswith("EX_") and met_mod1 in rxn.metabolites and len(rxn.metabolites) == 1
        ]
        if len(exchange_rxns) != 1:
            raise ValueError("There should be only one exchange reaction in model1 associated to"
                             " metabolite  " + met_mod1.id + " while there are "+ str(len(exchange_rxns)))
        ex_rxn = exchange_rxns[0]
        # Take bounds from the matching reaction.
        bounds = (ex_rxn.lower_bound, ex_rxn.upper_bound)
        # Remove exchange reaction
        cobra_model1.remove_reactions([ex_rxn])
        return bounds

    def create_exchange_reaction(self, copy_met_mod2: cobra.Metabolite, bounds):
        """
        Creates a new exchange reaction for the copy of met_mod2.
        Its ID is "EX_" + <processed met_mod2 id> + "_e".
        If bounds is provided as a tuple (lower_bound, upper_bound), those values are applied;
        otherwise, default bounds are used.
        """
        cobra_model1 = self.memo_model1.cobra_model
        new_ex_id = "EX_" + handle_metabolites_prefix_suffix(copy_met_mod2.id) + "_e"
        new_ex_rxn = cobra.Reaction(new_ex_id)
        new_ex_rxn.lower_bound, new_ex_rxn.upper_bound = bounds
        new_ex_rxn.add_metabolites({copy_met_mod2: -1})
        cobra_model1.add_reactions([new_ex_rxn])
        return new_ex_rxn

    def update_namespace_translation(self, met_mod2: cobra.Metabolite, met_mod1_translation_id: str) -> None:
        """
        Given a metabolite from model2 (met_mod2) and the expected translation ID (met_mod1_translation_id)
        for model1, this function:
          1. Finds the matching external metabolite in model1.
          2. Creates a translation reaction (TR_) to produce a copy of met_mod2.
          3. Deletes the corresponding exchange reaction for the original metabolite, capturing its bounds.
          4. Creates a new exchange reaction (EX_) for the new copy, carrying over the bounds.
        """
        # Step 1: Find the external metabolite in model1.
        # Fast lookup: Returns the external (compartment "e") metabolite in model1 whose processed ID matches the
        # given translation_id
        met_mod1 = self.external_metabolite_map.get(met_mod1_translation_id)
        if met_mod1 is None:
            msg = f"No external metabolite found in model1 for translation id '{met_mod1_translation_id}'"
            if self.fail_on_missing_metabolite:
                raise ValueError(msg)
            else:
                print(msg)
            return  # Skip creating the reaction for this metabolite.

        copy_met_mod2 = met_mod2.copy()

        # Step 2: Create the translation reaction.
        self.create_translation_reaction(met_mod1, copy_met_mod2)

        # Step 3: Delete the existing exchange reaction for met_mod1 and capture bounds.
        bounds = self.delete_exchange_reaction(met_mod1)

        # Step 4: Create a new exchange reaction for the copy of met_mod2 using the recovered bounds.
        self.create_exchange_reaction(copy_met_mod2, bounds)

        # Step 5: remove

    def translate_namespace(self):
        """
        Iterates over all external metabolites from model2 that are translatable
        (i.e. their processed IDs appear in the matches table's target_namespace) and updates model1's namespace.
        For each metabolite:
          - Retrieves the corresponding translation ID (source_namespace) from the matches table.
          - Calls update_namespace_translation to perform the translation reaction creation,
            delete the old exchange reaction, and create a new exchange reaction with carried-over bounds.
        """
        translatable_model2_metabs = [
            met for met in self.memo_model2.cobra_model.metabolites
            if met.compartment == "e" and handle_metabolites_prefix_suffix(met.id) in self.matches["target_namespace"].values
        ]
        for met_mod2 in translatable_model2_metabs:
            matching_row = self.matches[
                self.matches["target_namespace"] == handle_metabolites_prefix_suffix(met_mod2.id)
            ].iloc[0]
            met_mod1_translation_id = matching_row["source_namespace"]
            self.update_namespace_translation(met_mod2, met_mod1_translation_id)
