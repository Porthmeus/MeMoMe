from __future__ import annotations

from typing import Iterable

import cobra
import pandas as pd
import re
import warnings
from copy import deepcopy

from src import MeMoMetabolite
from src.MeMoModel import MeMoModel
from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix
from src.removeDuplicateMetabolites import removeDuplicateMetabolites





class ModelMerger:
    """
    merge two metabolic models by translating metabolites from one namespace to another
    using a translation compartment, and integrating shared metabolites.
    """
    def __init__(
        self,
        meMoModel1: MeMoModel,
        meMoModel2: MeMoModel,
        matches: pd.DataFrame,
        translation_compartment:str = "t"
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
        self.meMoModel1 = meMoModel1
        self.meMoModel2 = meMoModel2
        self.model1: cobra.Model = self.meMoModel1.cobra_model.copy()
        self.model2: cobra.Model = self.meMoModel2.cobra_model.copy()
        # remove duplicates and save cobra models on model1 and model2 - this
        # should only be done if the models have not been merged yet,
        # otherwise, we will remove the tranlation reactions
        # self.preprocess_models()
        self.matches = matches.copy()
        self.translation_compartment = translation_compartment
        if self.translation_compartment != "t":
            warnings.warn("Translation compartment uses " + self.translation_compartment + " instead of the standard 't' - be aware that this could lead to problems if further mergings should take place. Use consistent translation compartment ids across all merging attempts.")
        self.translated = False
        self.merged_model = None
    
    def copy_rxn(self, rxn):
        rxn_new = cobra.Reaction()
        rxn_new.annotation = deepcopy(rxn.annotation)
        rxn_new.bounds = deepcopy(rxn.bounds)
        #rxn_new.compartments = deepcopy(rxn.compartments)
        #rxn_new.gene_name_reaction_rule = rxn.gene_name_reaction_rule
        rxn_new.gene_reaction_rule = rxn.gene_reaction_rule
        #rxn_new.genes = deepcopy(rxn.genes)
        rxn_new.id = rxn.id
        #rxn_new.gpr = deepcopy(rxn.gpr)
        rxn_new.name = rxn.name
        rxn_new.notes = deepcopy(rxn.notes)
        #rxn_new.objective_coefficient = rxn.objective_coefficient
        rxn_new.reversibility = rxn.reversibility
        rxn_new.subsystem = rxn.subsystem
        rxn_new.add_metabolites(deepcopy(rxn.metabolites))
        return(rxn_new)

    def fetch_water_metabolites(self) -> [cobra.Metabolite, cobra.Metabolite]:
        # return a list of exactly two external water molecules for the two models to be merged
        ext_comp1 = cobra.medium.find_external_compartment(self.model1)
        ext_comp2 = cobra.medium.find_external_compartment(self.model2)
        ret = [None,None]
        for met1 in self.model1.metabolites:
            if (((met1.id.startswith("h2o|cpd00001|C00001|WATER")) or
                (met1.formula.lower() == "h2o") or 
                (met1.name.lower() == "water")) and
                (met1.compartment == ext_comp1)):
                    ret[0] = met1
                    break

        for met2 in self.model2.metabolites:
            if (((met2.id.startswith("h2o|cpd00001|C00001|WATER")) or
                (met2.formula.lower() == "h2o") or 
                (met2.name.lower() == "water")) and
                (met2.compartment == ext_comp2)):
                    ret[1] = met2
                    break

        if any([x == None for x in ret]):
            raise ValueError("I cannot find any water molecule for one of the models")
        return(ret)




    def add_translation_compartment(self, model:cb.Model):
        ''' Adds a new translation compartment in between the cell compartments and the exchange compartments. The metabolites will be translated via translation reactions which are preceded by an "TR_" prefix. This function will just add the layer but no translation yet '''
        # the order in which the metabolites and exchange reactions have to be
        # changed is a little awkward but otherwise we would need to create an
        # entirely new model, because the exchange reaction information is only
        # retained properly with this weird order.
        # 1. We create a copy of the exchange metabolite which will become the new exchange metabolite
        # 2. The original exchange metabolite becomes the translation metabolite (same namespace as the original model)
        # 3. We remove the translation metabolite (that was formerly the original exchange metabolite) from the exchange reaction
        # 4. We create a translation reaction, which transports the metabolite from the translation compartment to the external compartment

        ex_comp = cobra.medium.find_external_compartment(model)
        new_exchanges = []
        currents_mets = [x.id for x in model.metabolites]
        for met_id in currents_mets:
            met = model.metabolites.get_by_id(met_id)
            if met.compartment == ex_comp:
                # create a new metabolite which will become late the new exchange metabolite (1.) 
                met_newEx = met.copy()
                # rename all current exchange metabolites to translation metabolites (2.)
                met.compartment = self.translation_compartment
                met.id = re.sub(r"(_|\(|\[|__91__|__40__){compartment}($|\)$|\]$|__92__$|__41__$)".format(compartment = ex_comp),
                   r"\1{compartment}\2".format(compartment = self.translation_compartment),
                   met.id)
                ex_rxn = [x for x in met.reactions if x.id.startswith("EX_")][0]

                # standardize ex_rxn ids (so all end with _e)
                met_newEx.id = re.sub(r"(_|\(|\[|__91__|__40__){compartment}($|\)$|\]$|__92__$|__41__$)".format(compartment = ex_comp),
                   r"_{compartment}".format(compartment = ex_comp),
                   met_newEx.id)

                # remove translation metabolite from reaction and add the "new" external metabolite to the reaction (3.)
                ex_rxn.add_metabolites({met:1.0, met_newEx:-1.0})

                # add translation reaction (4.)
                tr_rxn = cobra.Reaction()
                tr_rxn.id = "TR_"+met.id
                tr_rxn.name = "Translation reaction for " + met.id
                tr_rxn.add_metabolites({met_newEx:1.0, met:-1.0}) 
                tr_rxn.lower_bound = -1000
                tr_rxn.upper_bound = 1000
                new_exchanges.append(tr_rxn)
        model.add_reactions(new_exchanges)
        model._compartments.update({self.translation_compartment:"translation"})
        return(model) 


    def translate(self) -> None:
        """Adds a translation compartment where metabolites are translated to the namespace of model1 """
        for mod in [self.model1, self.model2]:
            # create a new translation compartment
            if not self.translation_compartment in mod.compartments.keys():
                mod = self.add_translation_compartment(mod)
            # somebody else already used t as compartment
            else:
                if not "translation" in mod.compartments.values():
                    raise ValueError(f"Found a {self.translation_compartment} compartment in {mod.id} which is not used for translation - fix your models or change the translation compartment id while merging")
        
        # rename the metabolites for of the matching table in model2 and adjust stoichiometry for polymers

        ex_comp1 = cobra.medium.find_external_compartment(self.model1)
        ex_comp2 = cobra.medium.find_external_compartment(self.model2)
        exchange_ids2 = [x.id for x in self.model2.exchanges]
        h2o_1,h2o_2 = self.fetch_water_metabolites() 
        for met in self.model2.metabolites:
            if met.compartment == ex_comp2:
                met_id_basic = handle_metabolites_prefix_suffix(met.id)
                if met_id_basic in list(self.matches["met_id2"]):
                    new_met_id = self.matches.loc[self.matches["met_id2"] == met_id_basic,"met_id1"].item()
                    # first change exchange reaction name
                    for ex_rxn in met.reactions:
                        if ex_rxn.id in exchange_ids2:
                            ex_rxn.id = "EX_" + new_met_id + "_" + ex_comp1
                    # second change the metabolite name
                    met.id = new_met_id + "_" + ex_comp1
                    # third change compartment names
                    met.compartment = ex_comp1

        # forth adjust the stoichiometry to adjust for differing
        # polymer lengths
        # requires a new loop, otherwise we could add water before it is translated
        for met in self.model2.metabolites:
            if met.compartment == self.translation_compartment:
                met_id_basic = handle_metabolites_prefix_suffix(met.id)
                if met_id_basic in list(self.matches["met_id2"]):
                    new_met_id = self.matches.loc[self.matches["met_id2"] == met_id_basic,"met_id1"].item()
                    c_ratio = self.matches.loc[self.matches["met_id2"] == met_id_basic, "carbon_ratio"].item()
                    if not pd.isna(c_ratio) and c_ratio != 1:
                        # find the translation reaction
                        for tr_rxn in met.reactions:
                            if tr_rxn.id.startswith("TR_") and self.translation_compartment in tr_rxn.compartments:
                                mets_dic = tr_rxn.metabolites
                                met1 = list(mets_dic.keys())[0]
                                # add the stoichiometry - if the met1 is larger
                                if c_ratio > 1:
                                    new_dic = {met1 : 1,
                                               h2o_1 : (c_ratio-1),
                                               met : -c_ratio}
                                # add the stoichiometry - if the met2 is larger
                                elif c_ratio < 1:
                                    new_dic = {met1 : (1/c_ratio),
                                               h2o_1: -((1/c_ratio)-1),
                                               met: -1}
                                else:
                                    raise Error("Something is really wrong with the polymer stoichiometry setting")
                                # update the reaction
                                tr_rxn.subtract_metabolites(tr_rxn.metabolites)
                                tr_rxn.add_metabolites(new_dic)

        # finally add compartment annotation to the model
        self.model2._compartments.update({ex_comp1:self.model1.compartments[ex_comp1]})
        
        # set translated to true
        self.translated = True
        
    def merge_models(self) -> None:
        """ Merge both models together"""
        if self.translated == False: 
            self.translate()
        
        # check whether prefixes exist
        prefixes1 = self.model1.notes.get("MeMoMe_prefixes")
        prefixes2 = self.model2.notes.get("MeMoMe_prefixes")

        # harmonize prefixes if any of the models is already a merged model, otherwise add prefixes
        if prefixes1 == None :
            prefixed_model1 = self.add_prefix(self.model1.copy(), "M1")
            n_pre1 = 1
        else:
            n_pre1=len(prefixes1.keys())
            new_prefixes1 = ["M"+str(i+1) for i in range(n_pre1)]
            prefixed_model1 = self.correct_prefix(self.model1.copy(),
                                                  translation = dict(zip(list(prefixes1.keys()), new_prefixes1)))

        if prefixes2 == None:
            prefixed_model2 = self.add_prefix(self.model2.copy(), "M" + str(n_pre1+1))
        else:
            n_pre2=len(prefixes2.keys())
            new_prefixes2 = ["M"+str(i+1) for i in range(n_pre1,n_pre2)]
            prefixed_model2 = self.correct_prefix(self.model2.copy(),
                                                  translation = dict(zip(list(prefixes2.keys()), new_prefixes2)))
            #raise NotImplemented("Model prefixing currently works only for non-merged models")
        self.merge_models_simple(prefixed_model1, prefixed_model2)
    
    def merge_models_simple(self, prefixed_model1 : cobra.Model, prefixed_model2 : cobra.Model) -> None:
        """ Takes two prefixed models and simply merges them into one model """
        obj_rxns = [x.id for x in prefixed_model2.reactions if x.objective_coefficient == 1]
    
        # merge
        unique_prefixed2_rxns_ids = list({x.id for x in prefixed_model2.reactions} - {x.id for x in prefixed_model1.reactions})
        unique_prefixed2_rxns = [prefixed_model2.reactions.get_by_id(x) for x in unique_prefixed2_rxns_ids]
        prefixed_model1.add_reactions(unique_prefixed2_rxns)
        self.merged_model = prefixed_model1
        
        # adjust boundaries
        common_prefixed_rxns_ids = list({x.id for x in prefixed_model2.reactions} & {x.id for x in prefixed_model1.reactions})
        for rxn_id in common_prefixed_rxns_ids:
            lwbnd = min(prefixed_model1.reactions.get_by_id(rxn_id).lower_bound,
                        prefixed_model2.reactions.get_by_id(rxn_id).lower_bound)
            upbnd = max(prefixed_model1.reactions.get_by_id(rxn_id).upper_bound,
                        prefixed_model2.reactions.get_by_id(rxn_id).upper_bound)
            self.merged_model.reactions.get_by_id(rxn_id).lower_bound = lwbnd
            self.merged_model.reactions.get_by_id(rxn_id).upper_bound = upbnd

        # add objective coefficients
        for obj_r in obj_rxns:
            r = self.merged_model.reactions.get_by_id(obj_r)
            r.objective_coefficient = 1
        
        # change model id
        self.merged_model.id = "MeMoMe_merged_model"
        # add information for the prefixes to the model
        prefix_dict = deepcopy(prefixed_model1.notes["MeMoMe_prefixes"])
        prefix_dict.update(deepcopy(prefixed_model2.notes["MeMoMe_prefixes"]))
        self.merged_model.notes["MeMoMe_prefixes"] = prefix_dict

    def add_prefix(self, model: cobra.Model, prefix: str) -> cobra.Model:
        ''' Adds the prefix to the reactions and metabolites of the model (rxn1 -> prefix_rxn1) '''
        exch_comp = cobra.medium.find_external_compartment(model)
        exchange_rxns = [x.id for x in model.exchanges]
        external_mets = [x.id for x in model.metabolites if x.compartment == exch_comp]
        for rxn in model.reactions:
            if rxn.id not in exchange_rxns:
                rxn.id = prefix+ "_" + rxn.id
        for met in model.metabolites:
            if met.id not in external_mets:
                met.id = prefix + "_" + met.id
        model.notes["MeMoMe_prefixes"] = {prefix:model.id}
        return(model)
    
    def correct_prefix(self, model: cobra.Model, translation :dict) -> cobra.Model:
        ''' A simple correction for the prefixes '''
        # we need two rounds of corrections, to prevent overwriting of prefixes
        # (imagine you want to make M1 -> M2 and M2 -> M3, if you start with
        # M1, everthing will become M2, because in the second iteration the M2
        # will be converted to M3). Hence we paste old and new prefixes
        # together in the first round
        for old_prefix, new_prefix in translation.items():
            if not old_prefix == new_prefix:
                # replace prefixes in reactions
                for rxn in model.reactions:
                    if rxn.id.startswith(old_prefix + "_"):
                        rxn.id = new_prefix + rxn.id
                for met in model.metabolites:
                    if met.id.startswith(old_prefix+ "_"):
                        met.id = new_prefix + met.id

        # now we do the second round and correct the merge of both to the new one
        for old_prefix, new_prefix in translation.items():
            if not old_prefix == new_prefix:
                old_prefix = new_prefix + old_prefix + "_" 
                new_prefix = new_prefix + "_"
                # replace prefixes in reactions
                for rxn in model.reactions:
                    if rxn.id.startswith(old_prefix):
                        rxn.id = new_prefix + rxn.id[len(old_prefix):]
                for met in model.metabolites:
                    if met.id.startswith(old_prefix):
                        met.id = new_prefix + met.id[len(old_prefix):]
        old_prefix_notes = model.notes["MeMoMe_prefixes"]
        model.notes["MeMoMe_prefixes"] = dict(zip(translation.values(),old_prefix_notes.values()))
        return(model)
    
    def split_merged_model(self) -> list[cobra.Model]:
        ''' Take a the merged model and split it again to individual models '''
        if self.merged_model != None:
            mmod = self.merged_model.copy()
            ex_rxns = mmod.exchanges
            # get the objective reaction
            obj_rxn = [x.id for x in mmod.reactions if x.objective_coefficient !=0]
            # go through the different models and create a new model
            split_models = []
            for prefix,mod_name in mmod.notes["MeMoMe_prefixes"].items():
                # get the prefixed reactions 
                mod_rxns = [self.copy_rxn(x) for x in mmod.reactions if x.id.startswith(prefix + "_")]
                new_mod = cobra.Model()
                new_mod.add_reactions(mod_rxns)
                # remove the prefix and add objective coefficient
                for rxn in new_mod.reactions:
                    if rxn.id in obj_rxn:
                        rxn.objective_coefficient = 1
                    rxn.id = rxn.id[len(prefix) + 1:]
                for met in new_mod.metabolites:
                    if met.id.startswith(prefix + "_"):
                        met.id = met.id[len(prefix) + 1:]
                new_mod.add_reactions(ex_rxns)
                # adjust the id
                new_mod.id = mod_name
                split_models.append(new_mod)
            return(split_models)
        else:
            self.merge_models()
            self.split_merged_model()

    def preprocess_models(self) -> None:
        """Remove duplicate metabolites/reactions from both input models"""
        meMoModel1, _ = removeDuplicateMetabolites(self.meMoModel1)
        self.model1 = self.meMoModel1.cobra_model.copy()
        meMoModel2, _ = removeDuplicateMetabolites(self.meMoModel2)
        self.model2 = self.meMoModel2.cobra_model.copy()

