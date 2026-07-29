from __future__ import annotations

from typing import Iterable

import cobra
import pandas as pd
import re
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
        self.model1: cobra.Model | None  = None
        self.model2: cobra.Model | None  = None
        self.matches = matches.copy()
        self.translation_compartment = translation_compartment
        if self.translation_compartment != "t":
            warnings.warn("Translation compartment uses " + self.translation_compartment + " instead of the standard 't' - be aware that this could lead to problems if further mergings should take place. Use consistent translation compartment ids across all merging attempts.")
        # remove duplicates and save cobra models on model1 and model2
        self.preprocess_models()
        self.translated = False
    
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

    def add_translation_compartment(self, model:cb.Model):
        ex_comp = cobra.medium.find_external_compartment(model)
        new_exchanges = []
        currents_mets = [x.id for x in model.metabolites]
        for met_id in currents_mets:
            met = model.metabolites.get_by_id(met_id)
            if met.compartment == ex_comp:
                met_newEx = met.copy()
                # rename all current exchange metabolites to translation metabolites
                met.compartment = self.translation_compartment
                met.id = re.sub(r"(_|\(|\[|__91__|__40__){compartment}($|\)$|\]$|__92__$|__41__$)".format(compartment = ex_comp),
                   r"\1{compartment}\2".format(compartment = self.translation_compartment),
                   met.id)
                ex_rxn = [x for x in met.reactions if x.id.startswith("EX_")][0]
                met_newEx.id = re.sub(r"(_|\(|\[|__91__|__40__){compartment}($|\)$|\]$|__92__$|__41__$)".format(compartment = ex_comp),
                   r"_{compartment}".format(compartment = ex_comp),
                   met_newEx.id)
                ex_rxn.add_metabolites({met:1.0, met_newEx:-1.0}) # remove tranlation metabolite from reaction and add the "new" external metabolite to the reaction
                # standardize ex_rxn ids
                #ex_rxn.id = "EX_" + met_newEx.id

                # add translation reaction
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
        
        # rename the metabolites for of the matching table in model2
        ex_comp1 = cobra.medium.find_external_compartment(self.model1)
        ex_comp2 = cobra.medium.find_external_compartment(self.model2)
        exchange_ids2 = [x.id for x in self.model2.exchanges]
        for met in self.model2.metabolites:
            if met.compartment == ex_comp2:
                met_id_basic = handle_metabolites_prefix_suffix(met.id)
                new_met_id = self.matches.loc[self.matches["met_id2"] == met_id_basic,"met_id1"].item()
                # first change exchange reaction name
                for ex_rxn in met.reactions:
                    if ex_rxn.id in exchange_ids2:
                        ex_rxn.id = "EX_" + new_met_id + "_" + ex_comp1
                # second change the metabolite name
                met.id = new_met_id + "_" + ex_comp1
                # third change compartment names
                met.compartment = ex_comp1
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
            prefixed_model2 = self.add_prefix(self.model2.copy(), "M2")
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
                    if met.id.startswith(old_prefix):
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
        model.notes = dict(zip(translation.values(),old_prefix_notes.values()))
        return(model)


    def set_objective_reaction(self, target_model, merged_model: cobra.Model) -> None:
        objective_rxns = [rxn for rxn in target_model.reactions if rxn.objective_coefficient]
        if objective_rxns:
            original_id = objective_rxns[0].id
            # TODOs this fails for exchange reactions and model2 (or modelX)
            merged_obj_id = f"model1_{original_id}"
            if merged_obj_id in merged_model.reactions:
                merged_model.objective = merged_model.reactions.get_by_id(merged_obj_id)


    def preprocess_models(self) -> None:
        """Remove duplicate metabolites/reactions from both input models"""
        meMoModel1, _ = removeDuplicateMetabolites(self.meMoModel1)
        self.model1 = self.meMoModel1.cobra_model.copy()
        meMoModel2, _ = removeDuplicateMetabolites(self.meMoModel2)
        self.model2 = self.meMoModel2.cobra_model.copy()


 #   def add_prefix_to_model_ids(self, model: cobra.Model, prefix: str) -> None:
 #       """
 #       Adds a prefix to all metabolite and reaction IDs in a model to prevent conflicts.

 #       Parameters:
 #           model (cobra.Model): The model whose IDs will be prefixed.
 #           prefix (str): The prefix to add to the IDs.
 #       """
 #       for met in model.metabolites:
 #           met.id = prefix + met.id

 #       special_prefixes = ("EX_", "DM_", "SK_", "sink_")
 #       for reaction in model.reactions:
 #           original_id = reaction.id
 #           special_prefix = next((sp for sp in special_prefixes if original_id.startswith(sp)), None)
 #           if special_prefix:
 #               #TODO check if this senseful
 #               reaction.id = special_prefix + prefix + original_id[len(special_prefix) :]
 #           else:
 #               reaction.id = prefix + original_id

 #       for gene in model.genes:
 #           gene.id = f"{prefix}{gene.id}"

 #   # DONE 08_06_2026 (whole function)
 #   @staticmethod
 #   def is_external(metabolite: cobra.Metabolite) -> bool:
 #       return metabolite.compartment in ("e", "e0")
 #   
 #   # DONE 08_06_2026 (whole function)
 #   def is_water_metabolite(self, metabolite: MeMoMetabolite) -> bool:
 #       return ((metabolite._inchi_string == "InChI=1S/H2O/h1H2") or
 #               (metabolite._formula.lower() == "h2o") or (any("water" in name.lower() for name in metabolite._names))) and (self.is_external(metabolite))

 #   # DONE 08_06_2026 (whole function)
 #   def fetch_water_metabolites(self, meMoModel1:MeMoModel, meMoModel2:MeMoModel, model: cobra.Model) -> tuple[cobra.Metabolite, cobra.Metabolite]:
 #       """
 #       Fetch water metabolites from the merged cobra model for both source and target namespaces.
 #       If water metabolites are found in the original models, they should have been prefixed and added to the merged model. 
 #       This function identifies them based on their InChI, formula, or name, and returns the corresponding metabolite
 #       objects from the merged model. If no water metabolite is found for a namespace, the corresponding return value 
 #       will be None.

 #       Parameters:
 #           model (cobra.Model): The merged model containing both namespaces.
 #       """
 #       source_water_met = None
 #       target_water_met = None

 #       for memometab in meMoModel1.metabolites:
 #           #TODO:incorporate the is_external check in the is_water_metabolite function, and remove the is_external check from this loop
 #           if self.is_water_metabolite(memometab):
 #               # TODO model1 is hardcoded again
 #               target_water_met_id = f"model1_{memometab.id}"
 #               candidate_target_water_met = model.metabolites.get_by_id(target_water_met_id)
 #               if self.is_external(candidate_target_water_met):
 #                  target_water_met = candidate_target_water_met

 #       for memometab in meMoModel2.metabolites:
 #           if self.is_water_metabolite(memometab):
 #               source_water_met_id = f"model2_{memometab.id}"
 #               source_water_met = model.metabolites.get_by_id(source_water_met_id)
 #       if target_water_met is None or source_water_met is None:
 #           raise ValueError("No water metabolite found in one or both models.")
 #       
 #       return target_water_met, source_water_met
 #   
 #   def add_translation_metabolite(self, target_met_id: str) -> cobra.Metabolite:
 #       """
 #       Add (or fetch) a translation-compartment metabolite and its exchange reaction.

 #       Parameters:
 #           target_namespace (str): The target namespace identifier (without the _{self.TRANSLATION_COMPARTMENT} suffix).
 #       """
 #       def ensure_exchange_name(
 #           reaction: cobra.Reaction, translation_met: cobra.Metabolite
 #       ) -> None:
 #           if reaction.name:
 #               return
 #           met_name = translation_met.name or translation_met.id
 #           reaction.name = f"Translation compartment {met_name} exchange"

 #       translation_id = f"{target_met_id}_{self.TRANSLATION_COMPARTMENT}"
 #       # DONE 08_06_2026 (from this line to line 185 "self.merged_model.add_metabolites([translation_met])"
 #       source_met = self._find_source_metabolite(target_met_id)
 #       if translation_id in self.merged_model.metabolites:
 #           translation_met = self.merged_model.metabolites.get_by_id(translation_id)
 #           if source_met is not None:
 #               if translation_met.name is None:
 #                   translation_met.name = source_met.name
 #               if translation_met.formula is None:
 #                   translation_met.formula = source_met.formula
 #               if translation_met.charge is None:
 #                   translation_met.charge = source_met.charge
 #               if not translation_met.annotation and source_met.annotation:
 #                   translation_met.annotation = dict(source_met.annotation)

 #       if source_met is None:
 #           raise Exception(f"No source metabolite found for target ID {target_met_id!r}")
 #       else:
 #           name = source_met.name
 #           formula = source_met.formula
 #           charge = source_met.charge
 #           annotation = dict(source_met.annotation)
 #           translation_met = cobra.Metabolite(
 #               id=translation_id,
 #               name=name,
 #               formula=formula,
 #               charge=charge,
 #               compartment=self.TRANSLATION_COMPARTMENT
 #           )
 #           translation_met.annotation = annotation
 #           self.merged_model.add_metabolites([translation_met])

 #       ex_id = f"EX_{translation_id}"
 #       if ex_id not in self.merged_model.reactions:
 #           ex_rxn = cobra.Reaction(ex_id)
 #           ex_rxn.lower_bound = -1000
 #           ex_rxn.upper_bound = 1000
 #           ex_rxn.add_metabolites({translation_met: -1})
 #           ensure_exchange_name(ex_rxn, translation_met)
 #           self.merged_model.add_reactions([ex_rxn])
 #       else:
 #           ex_rxn = self.merged_model.reactions.get_by_id(ex_id)
 #           if translation_met not in ex_rxn.metabolites:
 #               ex_rxn.add_metabolites({translation_met: -1})
 #           ensure_exchange_name(ex_rxn, translation_met)
 #       return translation_met

 #   def _prefix_tr_reaction(self, reaction: cobra.Reaction, model_prefix: str) -> None:
 #       if model_prefix not in ("model1_", "model2_"):
 #           return
 #       if reaction.id.startswith(("TR_model1_", "TR_model2_")):
 #           return
 #       base_id = reaction.id
 #       if base_id.startswith("TR_"):
 #           base_id = base_id[len("TR_") :]
 #       if base_id.startswith(("model1_", "model2_")):
 #           base_id = base_id.split("_", 1)[1]
 #       model_tag = model_prefix.rstrip("_")
 #       reaction.id = f"TR_{model_tag}_{base_id}"

 #   # DONE 08_06_2026 (whole function, renamed target and added comments)
 #   def _find_source_metabolite(self, target_met_id: str) -> cobra.Metabolite | None:
 #       # Look for a metabolite object in the merged model that matches the target metabolite id, 
 #       # ignoring prefixes and suffixes, and prefer those in the extracellular compartment if 
 #       # multiple matching metabolites are found. Additionally, prefer a metabolite from the source model if possible. 
 #       # If not, return any match that is found, or None if no matches are found.
 #       def pick_from_prefix(prefix: str) -> cobra.Metabolite | None:
 #           fallback = None
 #           for met in self.merged_model.metabolites:
 #               if not met.id.startswith(prefix):
 #                   continue
 #               base_id = handle_metabolites_prefix_suffix(met.id[len(prefix) :])
 #               if base_id != target_met_id:
 #                   continue
 #               if self.is_external(met):
 #                   return met
 #               if fallback is None:
 #                   fallback = met
 #           return fallback

 #       return pick_from_prefix("model1_") or pick_from_prefix("model2_")
 #   
 #   def _build_exchange_map(
 #       self, model_prefix: str
 #   ) -> dict[str, list[tuple[cobra.Reaction, cobra.Metabolite]]]:
 #       exchange_map: dict[str, list[tuple[cobra.Reaction, cobra.Metabolite]]] = {}
 #       ex_prefix = f"EX_{model_prefix}"
 #       for rxn in self.merged_model.reactions:
 #           if not rxn.id.startswith(ex_prefix):
 #               continue
 #           for met in rxn.metabolites:
 #               if met.id.startswith(model_prefix):
 #                   base_id = handle_metabolites_prefix_suffix(met.id[len(model_prefix) :])
 #               else:
 #                   base_id = handle_metabolites_prefix_suffix(met.id)
 #               if base_id is None:
 #                   continue
 #               exchange_map.setdefault(base_id, []).append((rxn, met))
 #       return exchange_map

 #   def _propagate_subsystems_to_notes(self) -> None:
 #       """
 #       Persist subsystem information through SBML round-trips.

 #       cobra does not reliably write ``reaction.subsystem`` back into SBML,
 #       but it preserves ``reaction.notes``. To keep subsystem data, copy it
 #       into notes when missing.
 #       """
 #       for reaction in self.merged_model.reactions:
 #           subsystem = getattr(reaction, "subsystem", None)
 #           if not subsystem:
 #               continue
 #           notes = reaction.notes or {}
 #           if "SUBSYSTEM" not in notes or not notes.get("SUBSYSTEM"):
 #               notes["SUBSYSTEM"] = subsystem
 #           reaction.notes = notes


 #   def translate_namespace(
 #       self,
 #       #score_thr: float = 0.8,
 #       score_type: str = "total_score",
 #       inchi_score_thr: float = 1.0,
 #       name_score_thr: float = 0.9,
 #       db_score_thr: float = 0.5
 #   ) -> cobra.Model:
 #       """Public entry point that executes the full namespace translation pipeline."""
 #       target_exchange_map = self._build_exchange_map("model1_")
 #       source_exchange_map = self._build_exchange_map("model2_")
 #       if score_type not in self.matches.columns:
 #           raise ValueError(f"{score_type!r} not present in matches table columns")
 #       matches = self.matches[
 #           self.matches["target_namespace"].isin(target_exchange_map)
 #           & self.matches["source_namespace"].isin(source_exchange_map)
 #       ]
 #       matches = matches.loc[
 #           (matches["inchi_score"] == inchi_score_thr)
 #           | ((matches["Name_score"] >= name_score_thr) & (matches["DB_score"] >= db_score_thr))
 #       ]
##        if "Name_score" in matches.columns:
##            matches = matches.loc[matches["Name_score"] >= name_score_thr]
##        else:
##            raise ValueError("Name_score not present in matches table columns")
 #       matches = (
 #           matches.sort_values(
 #               by=[score_type, "target_namespace", "source_namespace"],
 #               ascending=[False, True, True],
 #           )
 #           .drop_duplicates(subset=["target_namespace"], keep="first")
 #       )

 #       # DONE 08_06_2026
 #       #fetch water metabolite for later use in polymerization or depolymerization of translation reactions, if needed
 #       target_water_met, source_water_met = self.fetch_water_metabolites(meMoModel1=self.meMoModel1, meMoModel2=self.meMoModel2, model=self.merged_model)

 #       for target,source in matches[['target_namespace', 'source_namespace']].itertuples(index=False):
 #           target_ex_rxn, target_met = target_exchange_map[target][0]
 #           # Add translation metabolite for each match, and creates its exchange reaction
 #           tr_met = self.add_translation_metabolite(target)
 #           source_ex_rxn, source_met = source_exchange_map[source][0]
 #           # Save original bounds
 #           source_lb = source_ex_rxn.lower_bound
 #           source_ub = source_ex_rxn.upper_bound
 #           target_lb = target_ex_rxn.lower_bound
 #           target_ub = target_ex_rxn.upper_bound
 #           # Rename exchange reactions replacing "EX_" with "TR_"
 #           source_ex_rxn.id = source_ex_rxn.id.replace("EX_", "TR_")
 #           target_ex_rxn.id = target_ex_rxn.id.replace("EX_", "TR_")
 #           self._prefix_tr_reaction(source_ex_rxn, "model2_")
 #           self._prefix_tr_reaction(target_ex_rxn, "model1_")
 #           # Set bounds to -1000 to 1000 to allow free flow
 #           source_ex_rxn.lower_bound = -1000
 #           source_ex_rxn.upper_bound = 1000
 #           target_ex_rxn.lower_bound = -1000
 #           target_ex_rxn.upper_bound = 1000
 #           # Modify source exchange reaction to export to translation compartment
 #           source_ex_rxn.add_metabolites({tr_met: 1})
 #           # Modify target exchange reaction to export to translation compartment
 #           target_ex_rxn.add_metabolites({tr_met: 1})
 #           # create ex_tr_met reaction to "import" from "nowhere" to target namespace.
 #           ex_tr_met_id = f"EX_{tr_met.id}"
 #           if ex_tr_met_id in self.merged_model.reactions:
 #               ex_tr_met_rxn = self.merged_model.reactions.get_by_id(ex_tr_met_id)
 #           else:
 #               ex_tr_met_rxn = cobra.Reaction(ex_tr_met_id)
 #               ex_tr_met_rxn.add_metabolites({tr_met: -1})
 #               self.merged_model.add_reactions([ex_tr_met_rxn])
 #           # TODO: convert the following bounds logic to a function
 #           # Set lower bound based on some function of source and target bounds (here we use sum as an example)
 #           ex_tr_met_rxn.lower_bound = source_lb + target_lb
 #           ex_tr_met_rxn.upper_bound = 1000

 #       # Convert any remaining exchange reactions (unmatched) to TR_ and add translation exchanges.
 #       for rxn in list(self.merged_model.reactions):
 #           if not rxn.id.startswith("EX_"):
 #               continue
 #           if any(met.compartment == self.TRANSLATION_COMPARTMENT for met in rxn.metabolites):
 #               continue

 #           met = None
 #           prefix = None
 #           for candidate in rxn.metabolites:
 #               if candidate.id.startswith("model1_"):
 #                   met = candidate
 #                   prefix = "model1_"
 #                   break
 #               if candidate.id.startswith("model2_"):
 #                   met = candidate
 #                   prefix = "model2_"
 #                   break
 #           if met is None:
 #               continue

 #           base_id = handle_metabolites_prefix_suffix(met.id[len(prefix) :])
 #           if base_id is None:
 #               continue

 #           tr_met = self.add_translation_metabolite(base_id)
 #           ex_lb = rxn.lower_bound
 #           ex_ub = rxn.upper_bound
 #           rxn.id = rxn.id.replace("EX_", "TR_", 1)
 #           self._prefix_tr_reaction(rxn, prefix)
 #           rxn.lower_bound = -1000
 #           rxn.upper_bound = 1000
 #           rxn.add_metabolites({tr_met: 1})

 #           ex_tr_met_id = f"EX_{tr_met.id}"
 #           if ex_tr_met_id in self.merged_model.reactions:
 #               ex_tr_met_rxn = self.merged_model.reactions.get_by_id(ex_tr_met_id)
 #               if tr_met not in ex_tr_met_rxn.metabolites:
 #                   ex_tr_met_rxn.add_metabolites({tr_met: -1})
 #           else:
 #               ex_tr_met_rxn = cobra.Reaction(ex_tr_met_id)
 #               ex_tr_met_rxn.add_metabolites({tr_met: -1})
 #               self.merged_model.add_reactions([ex_tr_met_rxn])
 #           ex_tr_met_rxn.lower_bound = ex_lb
 #           ex_tr_met_rxn.upper_bound = ex_ub

 #       # Ensure subsystem information survives SBML write/read cycles.
 #       self._propagate_subsystems_to_notes()
 #       self.set_objective_reaction(self.model1, self.merged_model)

 #       # Final cleanup: if any merged metabolite is missing a name, set it to the ID for better readability.
 #       for met in self.merged_model.metabolites:
 #           if met.name is None:
 #               met.name = met.id

 #       return self.merged_model
