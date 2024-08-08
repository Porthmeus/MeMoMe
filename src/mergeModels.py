import cobra
import re
import sys
import pandas as pd
from src.MeMoModel import *
sys.path.append('../')
from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix


class ModelMerger:
    """ The class responsible for managing the setting parameters for the merging procedure, including the list of 
    matching metabolites provided from the matching steps, as well as user settings like the definition of the new 
    common diet"""

    def __init__(self,
                 model1: MeMoModel,
                 model2: MeMoModel,
                 matches: pd.DataFrame
                 ):
        self.model1: MeMoModel = model1
        self.model2: MeMoModel = model2

        self.matches = matches

        self.model1 = self.step1(model1, "M1")
        self.model2 = self.step1(model2, "M2")

    def apply_translation(self, debug=False):
        if debug:
            reliable_matches = self.matches.loc[self.matches["Name_score"] == 1.0]
            print(reliable_matches)
            # print()
            # print(self.model1.cobra_model.compartments)
            # print()
            # print(self.model2.cobra_model.compartments)
            # print()

            # TODO: rename exchange metabolites and reactions of model2 to match the namespace of model1 (define step2 and call it)
            self.model2 = self.step2(model2, reliable_matches)
        else:
            raise NotImplemented()

        return

    def step1(self, me_mo_model: MeMoModel, prefix: str) -> MeMoModel:
        if "_" in prefix:
            raise ValueError("underscores are not allowed in prefix: " + prefix)
        # ensure that the model doesn't already have a common compartment
        if (len([reac.id for reac in me_mo_model.cobra_model.reactions if reac.id.startswith("EX_COMMON_")]) +
                len([met.id for met in me_mo_model.cobra_model.metabolites if met.id.startswith("COMMON_")]) == 0):
            # add prefix to reactions, metabolites and groups
            me_mo_model.cobra_model = self.adapt_ids(me_mo_model.cobra_model, prefix)
            # creates the common external metabolites and reactions and sets their upper and lower bounds
            me_mo_model.cobra_model = self.add_common(me_mo_model.cobra_model)

        else:
            raise NotImplemented()

        return me_mo_model

    def adapt_ids(self, model: cobra.Model, prefix: str) -> cobra.Model:
        """

        Parameters:
        model (cobra.Model): The COBRApy model object.
        prefix (str): The string prefix to prepend to all IDs.

        Returns:
        cobra.Model: The modified model with updated IDs.
        """
        ## Create a copy of the model to modify and return
        #model_copy = model.copy()

        # create a new compartment for the internal exchange
        model.compartments["i"] = "internal_exchange"

        # Update metabolite IDs adding prefix and modifying compartment
        for metabolite in model.metabolites:
            metabolite.id = prefix + "_" + metabolite.id
            if metabolite.compartment == "e":
                if metabolite.id.endswith("_e"):
                    metabolite.id = metabolite.id[:-2] + "_i"
                    metabolite.compartment = "i"
                # I change the compartment notation so from here on I can assume it is of this form
                elif metabolite.id.endswith("(e)"):
                    metabolite.id = metabolite.id[:-3] + "_i"
                elif metabolite.id.endswith("[e]"):
                    metabolite.id = metabolite.id[:-3] + "_i"
                else:
                    raise ValueError("The compartment specified in the metabolite id suffix does not correspond to the "
                                     "one specified in the compartment field for metabolite: " + metabolite.id)

        # replace the EX_ prefix with IEX and replace the compartment suffix
        for exch in model.exchanges:
            if exch.id.startswith("EX_"):
                # Replace 'EX_' with the new prefix
                exch.id = "IEX_" + prefix + "_" + exch.id[3:]
                if exch.id.endswith("_e"):
                    exch.id = exch.id[:-2] + "_i"
                # if the compartment's notation is not _e, I change it so from here on I can assume it is of this form
                elif exch.id.endswith("(e)"):
                    exch.id = exch.id[:-3] + "_i"
                elif exch.id.endswith("[e]"):
                    exch.id = exch.id[:-3] + "_i"
            else:
                raise ValueError("The exchange reaction's id: " + exch.id + "should start with the prefix 'EX_'")

        # Update gene IDs
        for gene in model.genes:
            gene.id = prefix + "_" + gene.id

        # Update group IDs
        if hasattr(model, 'groups'):
            for group in model.groups:
                group.id = prefix + "_" + group.id

        # Return the modified model
        return model

    def add_common(self, model: cobra.Model) -> cobra.Model:
        for r in [reac for reac in model.reactions if reac.id.startswith("IEX_")]:
            # Check if the reaction is an exchange reaction (only one metabolite)
            if len(r.metabolites) != 1:
                raise ValueError(f"The reaction {r.id} is not a valid exchange reaction.")

            # Get the original metabolite
            original_metabolite = list(r.metabolites.keys())[0]

            # Create a new metabolite with the same properties but different compartment
            # replace prefix and substitute the suffix to indicate the exchange compartment
            new_id = "COMMON"
            for el in original_metabolite.id.split('_')[1:-1]:
                new_id = new_id + "_" + el
            new_id = new_id + "_e"
            new_metabolite = original_metabolite.copy()
            new_metabolite.id = new_id
            new_metabolite.compartment = "e"


            # Add the new metabolite to the model
            model.add_metabolites([new_metabolite])
            # modify the IEX exchange reaction's stoichiometry in order to allow export to the COMMON compartment
            r.add_metabolites({new_metabolite: 1.0})
            # save lower and upper bound and make them unconstrained
            exch_lb = r.lower_bound
            exch_ub = r.upper_bound
            r.lower_bound = -1000
            r.upper_bound = 1000

        # Creates an exchange reaction for the new metabolite and adds it to the model
        model.add_boundary(new_metabolite, type="exchange", lb=exch_lb, ub=exch_ub)

        return model

    def step2(self, me_mo_model: MeMoModel, reliable_matches):
        for ex in me_mo_model.cobra_model.exchanges:
            #modify the metabolite id
            #print(ex)
            whole_met_id = list(ex.metabolites.keys())[0].id
            #print(whole_met_id)
            pattern = r"(?<=COMMON_).*?(?=_e)"
            met_id = re.search(pattern, whole_met_id)
            met_id = reliable_matches.loc[reliable_matches["met_id2"] == met_id]["met_id1"]
            new_met_id = "COMMON_" + met_id + "_e"
            list(ex.metabolites.keys())[0].id = new_met_id # update id with the translated metabolite id

            #TODO:modify the reaction id

        return


if __name__ == '__main__':
    tinyModel1 = "../tests/dat/my_model_with_annotations.xml"
    tinyModel2 = "../tests/dat/contiguous_pathway_model.xml"

    model1 = MeMoModel.fromPath(Path(tinyModel1))
    model2 = MeMoModel.fromPath(Path(tinyModel2))
    matches = model1.match(model2)
    print(matches)
    mo_me = ModelMerger(model1, model2, matches)
    mo_me.apply_translation(debug=True)
