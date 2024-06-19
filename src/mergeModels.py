import cobra
import pandas as pd
from src.MeMoModel import *


class ModelMerger:
    """ The class responsible for managing the setting parameters for the merging procedure, including the list of 
    matching metabolites provided from the matching steps, as well as user settings like the definition of the new 
    common diet"""

    def __init__(self,
                 model1: MeMoModel,
                 model2: MeMoModel,
                 matches: pd.DataFrame):
        self.model1: MeMoModel = model1
        self.model2: MeMoModel = model2
        self.matches: pd.DataFrame = matches

    def merge_models(self, debug=False):
        if debug:
            reliable_matches = self.matches.loc[self.matches["Name_score"] == 1.0]
            print(reliable_matches)
            print()
            print(self.model1.cobra_model.compartments)
            print()
            print(self.model2.cobra_model.compartments)
            print()
            self.model1 = self.step1(self.model1, prefix="M1_")
            # create an empty model
            merged_model = MeMoModel(cobra_model=cobra.Model(id_or_model="merged_model", name="merged_model"))

        else:
            raise NotImplemented()

        return merged_model

    def adapt_ids(model: cobra.Model, prefix: str) -> cobra.Model:
        """

        Parameters:
        model (cobra.Model): The COBRApy model object.
        prefix (str): The string prefix to prepend to all IDs.

        Returns:
        cobra.Model: The modified model with updated IDs.
        """
        # Create a copy of the model to modify and return
        model_copy = model.copy()

        # rename the external compartment to internal
        model_copy.compartments.pop("e")
        model_copy.compartments["i"] = "internal_exchange"

        # Update metabolite IDs adding prefix and modifying compartment
        for metabolite in model_copy.metabolites:
            metabolite.id = prefix + metabolite.id
            if metabolite.id.endswith("_e"):
                if metabolite.compartment == "e":
                    metabolite.id = metabolite.id[:-2] + "_i"
                else:
                    raise ValueError("The compartment specified in the metabolite id suffix does not correspond to the "
                                     "one specified in the compartment field for metabolite: " + metabolite.id)

        # replace the EX_ prefix with IEX and replace the compartment suffix
        for reaction in model_copy.reactions:
            reaction.id = reaction.id.replace("EX_")
            reaction.id = prefix + reaction.id

        # Update gene IDs
        for gene in model_copy.genes:
            gene.id = prefix + gene.id

        # Update group IDs
        if hasattr(model_copy, 'groups'):
            for group in model_copy.groups:
                group.id = prefix + group.id

        # Return the modified model
        return model_copy

    def step1(self, model: MeMoModel, prefix: str) -> MeMoModel:
        # ensure that the model doesn't already have a common compartment
        if (len([reac.id for reac in model.cobra_model.reactions if reac.id.startswith("EX_COMMON_")]) +
                len([met.id for met in model.cobra_model.metabolites if met.id.startswith("COMMON_")]) == 0):
            # add prefix to reactions, metabolites and groups
            model.cobra_model = self.adapt_ids(model.cobra_model, prefix)
            # add the feeding compartment

            # add reaction to move from the external to the feeding compartment
            # "move" the lower bounds from the external to the feeding compartments

        else:
            raise NotImplemented()

        return model

if __name__ == '__main__':
    tinyModel1 = "../tests/dat/my_model_with_annotations.xml"
    tinyModel2 = "../tests/dat/contiguous_pathway_model.xml"

    model1 = MeMoModel.fromPath(Path(tinyModel1))
    model2 = MeMoModel.fromPath(Path(tinyModel2))
    matches = model1.match(model2)
    mo_me = ModelMerger(model1, model2, matches)
    new_model = mo_me.merge_models(debug=True)
