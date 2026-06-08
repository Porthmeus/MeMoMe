from __future__ import annotations

from typing import Iterable

import cobra
import pandas as pd

from src import MeMoMetabolite
from src.MeMoModel import MeMoModel
from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix
from src.removeDuplicateMetabolites import removeDuplicateMetabolites




class ModelMerger:
    """
    merge two metabolic models by translating metabolites from one namespace to another
    using a translation compartment, and integrating shared metabolites.
    """
    TRANSLATION_COMPARTMENT = "t"
    def __init__(
        self,
        meMoModel1: MeMoModel,
        meMoModel2: MeMoModel,
        matches: pd.DataFrame,
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
        self.matches = matches.copy()
        # Rename columns for clarity
        self.matches = self.matches.rename(columns={'met_id2': 'source_namespace'})
        self.matches = self.matches.rename(columns={'met_id1': 'target_namespace'})
        # Ensure that the matches DataFrame contains necessary columns
        required_columns = {'source_namespace', 'target_namespace'}
        if not required_columns.issubset(self.matches.columns):
            raise ValueError(f"The matches DataFrame must contain columns: {required_columns}")
        # Preprocess models to handle duplicates
        self.preprocess_models()
        self.add_prefix_to_model_ids(self.model1, "model1_")
        self.add_prefix_to_model_ids(self.model2, "model2_")
        # No need to update matches DataFrame with prefixes since we'll remove them during translation
        # Initialize merged model
        self.merged_model = cobra.Model("merged_model")
        # Add metabolites and reactions from both models, cobra.Model does not expose add_genes; genes referenced in GPRs are already added via add_reactions
        self.merged_model.add_metabolites(self.model1.metabolites)
        self.merged_model.add_reactions(self.model1.reactions)
        self.merged_model.add_metabolites(self.model2.metabolites)
        self.merged_model.add_reactions(self.model2.reactions)


    def set_objective_reaction(self, target_model, merged_model: cobra.Model) -> None:
        objective_rxns = [rxn for rxn in target_model.reactions if rxn.objective_coefficient]
        if objective_rxns:
            original_id = objective_rxns[0].id
            merged_obj_id = f"model1_{original_id}"
            if merged_obj_id in merged_model.reactions:
                merged_model.objective = merged_model.reactions.get_by_id(merged_obj_id)


    def preprocess_models(self) -> None:
        """Remove duplicate metabolites/reactions from both input models"""
        #meMoModel1 = self.meMoModel1#MeMoModel(cobra_model=self.model1)
        #meMoModel2 = self.meMoModel2#MeMoModel(cobra_model=self.model2)
        meMoModel1, _ = removeDuplicateMetabolites(self.meMoModel1)
        self.model1 = self.meMoModel1.cobra_model.copy()
        meMoModel2, _ = removeDuplicateMetabolites(self.meMoModel2)
        self.model2 = self.meMoModel2.cobra_model.copy()


    def add_prefix_to_model_ids(self, model: cobra.Model, prefix: str) -> None:
        """
        Adds a prefix to all metabolite and reaction IDs in a model to prevent conflicts.

        Parameters:
            model (cobra.Model): The model whose IDs will be prefixed.
            prefix (str): The prefix to add to the IDs.
        """
        for met in model.metabolites:
            met.id = prefix + met.id

        special_prefixes = ("EX_", "DM_", "SK_", "sink_")
        for reaction in model.reactions:
            original_id = reaction.id
            special_prefix = next((sp for sp in special_prefixes if original_id.startswith(sp)), None)
            if special_prefix:
                reaction.id = special_prefix + prefix + original_id[len(special_prefix) :]
            else:
                reaction.id = prefix + original_id

        for gene in model.genes:
            gene.id = f"{prefix}{gene.id}"

    # DONE 08_06_2026 (whole function)
    @staticmethod
    def is_external(metabolite: cobra.Metabolite) -> bool:
        return metabolite.compartment in ("e", "e0")
    
    # DONE 08_06_2026 (whole function)
    def is_water_metabolite(self, metabolite: MeMoMetabolite, model: cobra.Model) -> bool:
        return ((metabolite._inchi == "InChI=1S/H2O/h1H2") or (metabolite._formula.lower() == "h2o")
                 or (any("water" in name.lower() for name in metabolite._names))) and (self.is_external(metabolite, model))

    # DONE 08_06_2026 (whole function)
    def fetch_water_metabolites(self, memomodel1:MeMoModel, memomodel2:MeMoModel, model: cobra.Model) -> tuple[cobra.Metabolite, cobra.Metabolite]:
        """
        Fetch water metabolites from the merged cobra model for both source and target namespaces.
        If water metabolites are found in the original models, they should have been prefixed and added to the merged model. 
        This function identifies them based on their InChI, formula, or name, and returns the corresponding metabolite
        objects from the merged model. If no water metabolite is found for a namespace, the corresponding return value 
        will be None.

        Parameters:
            model (cobra.Model): The merged model containing both namespaces.
        """
        source_water_met = None
        target_water_met = None

        for memometab in memomodel1.metabolites:
            #TODO:incorporate the is_external check in the is_water_metabolite function, and remove the is_external check from this loop
            if self.is_water_metabolite(memometab):
                target_water_met_id = f"model1_{memometab.id}"
                candidate_target_water_met = model.metabolites.get_by_id(target_water_met_id)
                if self.is_external(candidate_target_water_met):
                   target_water_met = candidate_target_water_met

        for memometab in memomodel2.metabolites:
            if self.is_water_metabolite(memometab):
                source_water_met_id = f"model2_{memometab.id}"
                source_water_met = model.metabolites.get_by_id(source_water_met_id)
        if target_water_met is None or source_water_met is None:
            raise ValueError("No water metabolite found in one or both models.")
        
        return target_water_met, source_water_met
    
    def add_translation_metabolite(self, target_met_id: str) -> cobra.Metabolite:
        """
        Add (or fetch) a translation-compartment metabolite and its exchange reaction.

        Parameters:
            target_namespace (str): The target namespace identifier (without the _{self.TRANSLATION_COMPARTMENT} suffix).
        """
        def ensure_exchange_name(
            reaction: cobra.Reaction, translation_met: cobra.Metabolite
        ) -> None:
            if reaction.name:
                return
            met_name = translation_met.name or translation_met.id
            reaction.name = f"Translation compartment {met_name} exchange"

        translation_id = f"{target_met_id}_{self.TRANSLATION_COMPARTMENT}"
        # DONE 08_06_2026 (from this line to line 185 "self.merged_model.add_metabolites([translation_met])"
        source_met = self._find_source_metabolite(target_met_id)
        if translation_id in self.merged_model.metabolites:
            translation_met = self.merged_model.metabolites.get_by_id(translation_id)
            if source_met is not None:
                if translation_met.name is None:
                    translation_met.name = source_met.name
                if translation_met.formula is None:
                    translation_met.formula = source_met.formula
                if translation_met.charge is None:
                    translation_met.charge = source_met.charge
                if not translation_met.annotation and source_met.annotation:
                    translation_met.annotation = dict(source_met.annotation)

        if source_met is None:
            raise Exception(f"No source metabolite found for target ID {target_met_id!r}")
        else:
            name = source_met.name
            formula = source_met.formula
            charge = source_met.charge
            annotation = dict(source_met.annotation)
            translation_met = cobra.Metabolite(
                id=translation_id,
                name=name,
                formula=formula,
                charge=charge,
                compartment=self.TRANSLATION_COMPARTMENT
            )
            translation_met.annotation = annotation
            self.merged_model.add_metabolites([translation_met])

        ex_id = f"EX_{translation_id}"
        if ex_id not in self.merged_model.reactions:
            ex_rxn = cobra.Reaction(ex_id)
            ex_rxn.lower_bound = -1000
            ex_rxn.upper_bound = 1000
            ex_rxn.add_metabolites({translation_met: -1})
            ensure_exchange_name(ex_rxn, translation_met)
            self.merged_model.add_reactions([ex_rxn])
        else:
            ex_rxn = self.merged_model.reactions.get_by_id(ex_id)
            if translation_met not in ex_rxn.metabolites:
                ex_rxn.add_metabolites({translation_met: -1})
            ensure_exchange_name(ex_rxn, translation_met)
        return translation_met

    def _prefix_tr_reaction(self, reaction: cobra.Reaction, model_prefix: str) -> None:
        if model_prefix not in ("model1_", "model2_"):
            return
        if reaction.id.startswith(("TR_model1_", "TR_model2_")):
            return
        base_id = reaction.id
        if base_id.startswith("TR_"):
            base_id = base_id[len("TR_") :]
        if base_id.startswith(("model1_", "model2_")):
            base_id = base_id.split("_", 1)[1]
        model_tag = model_prefix.rstrip("_")
        reaction.id = f"TR_{model_tag}_{base_id}"

    # DONE 08_06_2026 (whole function, renamed target and added comments)
    def _find_source_metabolite(self, target_met_id: str) -> cobra.Metabolite | None:
        # Look for a metabolite object in the merged model that matches the target metabolite id, 
        # ignoring prefixes and suffixes, and prefer those in the extracellular compartment if 
        # multiple matching metabolites are found. Additionally, prefer a metabolite from the source model if possible. 
        # If not, return any match that is found, or None if no matches are found.
        def pick_from_prefix(prefix: str) -> cobra.Metabolite | None:
            fallback = None
            for met in self.merged_model.metabolites:
                if not met.id.startswith(prefix):
                    continue
                base_id = handle_metabolites_prefix_suffix(met.id[len(prefix) :])
                if base_id != target_met_id:
                    continue
                if self.is_external(met):
                    return met
                if fallback is None:
                    fallback = met
            return fallback

        return pick_from_prefix("model1_") or pick_from_prefix("model2_")
    
    def _build_exchange_map(
        self, model_prefix: str
    ) -> dict[str, list[tuple[cobra.Reaction, cobra.Metabolite]]]:
        exchange_map: dict[str, list[tuple[cobra.Reaction, cobra.Metabolite]]] = {}
        ex_prefix = f"EX_{model_prefix}"
        for rxn in self.merged_model.reactions:
            if not rxn.id.startswith(ex_prefix):
                continue
            for met in rxn.metabolites:
                if met.id.startswith(model_prefix):
                    base_id = handle_metabolites_prefix_suffix(met.id[len(model_prefix) :])
                else:
                    base_id = handle_metabolites_prefix_suffix(met.id)
                if base_id is None:
                    continue
                exchange_map.setdefault(base_id, []).append((rxn, met))
        return exchange_map

    def _propagate_subsystems_to_notes(self) -> None:
        """
        Persist subsystem information through SBML round-trips.

        cobra does not reliably write ``reaction.subsystem`` back into SBML,
        but it preserves ``reaction.notes``. To keep subsystem data, copy it
        into notes when missing.
        """
        for reaction in self.merged_model.reactions:
            subsystem = getattr(reaction, "subsystem", None)
            if not subsystem:
                continue
            notes = reaction.notes or {}
            if "SUBSYSTEM" not in notes or not notes.get("SUBSYSTEM"):
                notes["SUBSYSTEM"] = subsystem
            reaction.notes = notes


    def translate_namespace(
        self,
        #score_thr: float = 0.8,
        score_type: str = "total_score",
        inchi_score_thr: float = 1.0,
        name_score_thr: float = 0.9,
        db_score_thr: float = 0.5
    ) -> cobra.Model:
        """Public entry point that executes the full namespace translation pipeline."""
        target_exchange_map = self._build_exchange_map("model1_")
        source_exchange_map = self._build_exchange_map("model2_")
        if score_type not in self.matches.columns:
            raise ValueError(f"{score_type!r} not present in matches table columns")
        matches = self.matches[
            self.matches["target_namespace"].isin(target_exchange_map)
            & self.matches["source_namespace"].isin(source_exchange_map)
        ]
        matches = matches.loc[
            (matches["inchi_score"] == inchi_score_thr)
            | ((matches["Name_score"] >= name_score_thr) & (matches["DB_score"] >= db_score_thr))
        ]
#        if "Name_score" in matches.columns:
#            matches = matches.loc[matches["Name_score"] >= name_score_thr]
#        else:
#            raise ValueError("Name_score not present in matches table columns")
        matches = (
            matches.sort_values(
                by=[score_type, "target_namespace", "source_namespace"],
                ascending=[False, True, True],
            )
            .drop_duplicates(subset=["target_namespace"], keep="first")
        )

        # DONE 08_06_2026
        #fetch water metabolite for later use in polymerization or depolymerization of translation reactions, if needed
        target_water_met, source_water_met = self.fetch_water_metabolites(memomodel1=self.memomodel1, memomodel2=self.memomodel2, model=self.merged_model)

        for target,source in matches[['target_namespace', 'source_namespace']].itertuples(index=False):
            target_ex_rxn, target_met = target_exchange_map[target][0]
            # Add translation metabolite for each match, and creates its exchange reaction
            tr_met = self.add_translation_metabolite(target)
            source_ex_rxn, source_met = source_exchange_map[source][0]
            # Save original bounds
            source_lb = source_ex_rxn.lower_bound
            source_ub = source_ex_rxn.upper_bound
            target_lb = target_ex_rxn.lower_bound
            target_ub = target_ex_rxn.upper_bound
            # Rename exchange reactions replacing "EX_" with "TR_"
            source_ex_rxn.id = source_ex_rxn.id.replace("EX_", "TR_")
            target_ex_rxn.id = target_ex_rxn.id.replace("EX_", "TR_")
            self._prefix_tr_reaction(source_ex_rxn, "model2_")
            self._prefix_tr_reaction(target_ex_rxn, "model1_")
            # Set bounds to -1000 to 1000 to allow free flow
            source_ex_rxn.lower_bound = -1000
            source_ex_rxn.upper_bound = 1000
            target_ex_rxn.lower_bound = -1000
            target_ex_rxn.upper_bound = 1000
            # Modify source exchange reaction to export to translation compartment
            source_ex_rxn.add_metabolites({tr_met: 1})
            # Modify target exchange reaction to export to translation compartment
            target_ex_rxn.add_metabolites({tr_met: 1})
            # create ex_tr_met reaction to "import" from "nowhere" to target namespace.
            ex_tr_met_id = f"EX_{tr_met.id}"
            if ex_tr_met_id in self.merged_model.reactions:
                ex_tr_met_rxn = self.merged_model.reactions.get_by_id(ex_tr_met_id)
            else:
                ex_tr_met_rxn = cobra.Reaction(ex_tr_met_id)
                ex_tr_met_rxn.add_metabolites({tr_met: -1})
                self.merged_model.add_reactions([ex_tr_met_rxn])
            # TODO: convert the following bounds logic to a function
            # Set lower bound based on some function of source and target bounds (here we use sum as an example)
            ex_tr_met_rxn.lower_bound = source_lb + target_lb
            ex_tr_met_rxn.upper_bound = 1000

        # Convert any remaining exchange reactions (unmatched) to TR_ and add translation exchanges.
        for rxn in list(self.merged_model.reactions):
            if not rxn.id.startswith("EX_"):
                continue
            if any(met.compartment == self.TRANSLATION_COMPARTMENT for met in rxn.metabolites):
                continue

            met = None
            prefix = None
            for candidate in rxn.metabolites:
                if candidate.id.startswith("model1_"):
                    met = candidate
                    prefix = "model1_"
                    break
                if candidate.id.startswith("model2_"):
                    met = candidate
                    prefix = "model2_"
                    break
            if met is None:
                continue

            base_id = handle_metabolites_prefix_suffix(met.id[len(prefix) :])
            if base_id is None:
                continue

            tr_met = self.add_translation_metabolite(base_id)
            ex_lb = rxn.lower_bound
            ex_ub = rxn.upper_bound
            rxn.id = rxn.id.replace("EX_", "TR_", 1)
            self._prefix_tr_reaction(rxn, prefix)
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000
            rxn.add_metabolites({tr_met: 1})

            ex_tr_met_id = f"EX_{tr_met.id}"
            if ex_tr_met_id in self.merged_model.reactions:
                ex_tr_met_rxn = self.merged_model.reactions.get_by_id(ex_tr_met_id)
                if tr_met not in ex_tr_met_rxn.metabolites:
                    ex_tr_met_rxn.add_metabolites({tr_met: -1})
            else:
                ex_tr_met_rxn = cobra.Reaction(ex_tr_met_id)
                ex_tr_met_rxn.add_metabolites({tr_met: -1})
                self.merged_model.add_reactions([ex_tr_met_rxn])
            ex_tr_met_rxn.lower_bound = ex_lb
            ex_tr_met_rxn.upper_bound = ex_ub

        # Ensure subsystem information survives SBML write/read cycles.
        self._propagate_subsystems_to_notes()
        self.set_objective_reaction(self.model1, self.merged_model)

        # Final cleanup: if any merged metabolite is missing a name, set it to the ID for better readability.
        for met in self.merged_model.metabolites:
            if met.name is None:
                met.name = met.id

        return self.merged_model
