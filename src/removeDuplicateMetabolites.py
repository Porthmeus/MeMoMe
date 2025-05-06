# Porthmeus
# 28.05.24

import re
import cobra as cb
import copy
import pandas as pd
from src.MeMoModel import MeMoModel


def printModelStats(model: cb.Model):
    print(f"{len(model.reactions)} reactions")
    print(f"{len(model.metabolites)} metabolites")
    print(f"{len(model.genes)} genes")


def findCommonReactions(
    met1: cb.Metabolite, met2: cb.Metabolite, reversible_is_same=True
) -> list[cb.Reaction]:
    """Takes two metabolites of the same model and returns reactions which are the same except for the two metabolites
    met1,met2 - the two cobra metabolite objects which are the categorized as same
    reverse_is_same - boolean, should reactions be reported if one of them is reversible but the other is not
    """

    # check if the metabolites originate from the same model
    if met1.model != met2.model:
        raise ValueError("Metabolites do not originate from the same model! Aborted.")

    # get the essential data
    reactions1 = [x for x in met1.reactions]
    reactions2 = [x for x in met2.reactions]
    reaction_rules1 = [x.reaction for x in reactions1]
    reaction_rules2 = [x.reaction for x in reactions2]
    metId1 = met1.id
    metId2 = met2.id

    # iterate through list of reactions
    same = []
    for i in range(len(reactions1)):
        for j in range(len(reactions2)):
            reac_rule1 = copy.deepcopy(reaction_rules1[i])
            reac_rule2 = copy.deepcopy(reaction_rules2[j])
            if reversible_is_same == True:
                if reac_rule1.find("<=>") != -1 or reac_rule2.find("<=>") != -1:
                    reac_rule1 = reac_rule1.replace("-->", "<=>").replace("<--", "<=>")
                    reac_rule2 = reac_rule2.replace("-->", "<=>").replace("<--", "<=>")
            # else:
            # print(reac_rule1)
            # print(reac_rule2)
            # TODO: check if compartments are the same for the reactions. Kind of difficult as several compartments are allowed for one reaction - not sure how to handle that for now
            if reac_rule1 == re.sub(metId2, metId1, reac_rule2):
                same.append((reactions1[i], reactions2[j]))
    return same


def mergeReactions(
    rxn1: cb.Reaction, rxn2: cb.Reaction, remove_orphans: bool = True, verbose=False
) -> None:
    """Takes a pair of reactions from the same model and merges them. This will result in merging the information into the first reaction and deletion of the second.

    remove_orphans = boolean, if True (default) - all orphanized metabolites and genes in the model will be deleted
    """

    # check if compartments are the same
    if len(rxn1.compartments) > 1:
        raise NotImplementedError(
            "Reaction "
            + rxn1.id
            + " is associated to more than one compartment, don't know how to handle that"
        )
    if len(rxn2.compartments) > 1:
        raise NotImplementedError(
            "Reaction "
            + rxn2.id
            + " is associated to more than one compartment, don't know how to handle that"
        )
    assert (
        rxn1.compartments == rxn2.compartments
    ), f"You try to merge reactions which take place in different compartments ({rxn1.id=}, {rxn2.id =})"

    # check if the gprs are the same, if not merge
    if rxn1.gene_reaction_rule != rxn2.gene_reaction_rule:
        rxn1.gene_reaction_rule = (
            "(" + rxn1.gene_reaction_rule + ") or (" + rxn2.gene_reaction_rule + ")"
        )

    # find the most liberal bounds
    rxn1.lower_bound = min([rxn1.lower_bound, rxn2.lower_bound])
    rxn1.upper_bound = max([rxn1.upper_bound, rxn2.upper_bound])

    # merge annotations
    # check if there is annotations at all
    if rxn1.annotation != None and rxn2.annotation:
        annotation = rxn1.annotation
        annotation2 = rxn2.annotation

        # convert all entries into lists
        for x in annotation.keys():
            if type(annotation[x]) != list:
                annotation[x] = [annotation[x]]
        for x in annotation2.keys():
            if type(annotation2[x]) != list:
                annotation2[x] = [annotation2[x]]

        # merge the annotation
        for x in annotation2.keys():
            if x in annotation.keys():
                annotation[x].extend(annotation2[x])
                # remove duplicates and sort
                new_annolst = list(set(annotation[x]))
                new_annolst.sort()
                self.annotations[x] = new_annolst
    elif rxn1.annotation != None:
        annotation = rxn1.annotation
    elif rxn2.annotation != None:
        annotation = rxn2.annotation
    else:
        annotation = dict()
    rxn1.annotation = annotation

    # all information of reaction 2 should now be included in reaction 1, thus we can safely remove the second reaction
    rxn2.remove_from_model(remove_orphans=remove_orphans)
    if verbose:
        print(f"Removed {rxn2.id=} from {rxn1.model.id=}")


def exchangeMetabolite(
    met1: cb.Metabolite, met2: cb.Metabolite, prune: bool = True
) -> pd.DataFrame:
    """Find the all reactions which use the second metabolite and exchange it with the first metabolite"""

    rxns = met2.reactions
    # keep the information which reactions are affected
    rxns_ids = [x.id for x in rxns]
    for rxn in rxns:
        coef = rxn.metabolites[met2]
        sub_dic = {met2: coef, met1: -1 * coef}
        rxn.subtract_metabolites(sub_dic)

    if prune == True:
        met2.remove_from_model()

    # cast the information into a pandas data frame
    n = len(rxns_ids)
    res = pd.DataFrame(
        {
            "new_rxn_id": rxns_ids,
            "old_rxn_id": rxns_ids,
            "new_met_id": [met1.id] * len(rxns_ids),
            "old_met_id": [met2.id] * len(rxns_ids),
            "what": ["met_exchange"] * len(rxns_ids),
        }
    )
    return res


def mergeAndRemove(
    met1: cb.Metabolite,
    met2: cb.Metabolite,
    prune: bool = True,
    reversible_is_same: bool = True,
) -> pd.DataFrame:
    """Clean the model of a duplicated metabolite. Every instance of met2 will be replaced by met1. Reactions which essentilly have the same"""
    # first find duplicated reactions and merge them
    rxns = findCommonReactions(met1, met2, reversible_is_same=reversible_is_same)
    # ceate a data frame to keep track of the changes
    res1 = pd.DataFrame(
        {
            "new_rxn_id": [x[0].id for x in rxns],
            "old_rxn_id": [x[1].id for x in rxns],
            "new_met_id": [met1.id] * len(rxns),
            "old_met_id": [met2.id] * len(rxns),
            "what": ["rxn_merge"] * len(rxns),
        }
    )
    for rxn_pair in rxns:
        mergeReactions(rxn_pair[0], rxn_pair[1])
    # exchange duplicated metabolite
    res2 = exchangeMetabolite(met1, met2, prune=prune)
    res = pd.concat([res1, res2])
    return res


def detectDuplicates(selfmatch: pd.DataFrame, check_charge: bool = True) -> list[tuple]:
    """Takes the result from selfmatch of a model and tries to detect all duplicated metabolites. Returns a list of tuples with the duplicated but not identical metabolites"""
    # some renaming
    matches = selfmatch

    # get the top score for each metabolite
    matches_max1 = pd.DataFrame(matches.groupby("met_id1")["total_score"].max())
    matches_max2 = pd.DataFrame(matches.groupby("met_id2")["total_score"].max())

    # prefilter the data - this has performance reasons, but I could not figure out how to do the actual filtering in a vectorized way
    # check if it is not a self match
    matches = matches.loc[matches["met_id1"] != matches["met_id2"]]
    # check whether the first metabolite has a match on the highest total score
    matches.index = matches["met_id1"]
    keep_max1 = matches.eq(matches_max1)["total_score"]
    # check if the second metabolite has reached the maximum score
    matches.index = matches["met_id2"]
    keep_max2 = matches.eq(matches_max2)["total_score"]
    keep = [
        x for x, y in enumerate(keep_max1) if y == True and keep_max2.iloc[x] == True
    ]
    met_id1 = keep_max1.iloc[keep].index
    met_id2 = keep_max2.iloc[keep].index
    matches = matches.loc[matches["met_id1"].isin(met_id1)]
    matches = matches.loc[matches["met_id2"].isin(met_id2)]

    # do the actual filtering iteratively
    matches = matches.reset_index(drop=True)  # reset index
    keep = []
    for i, row in matches.iterrows():
        if all(
            [
                # check whether the first metabolite has a match on the highest total score
                matches_max1["total_score"][row["met_id1"]] == row["total_score"],
                # check if the second metabolite has reached the maximum score
                matches_max2["total_score"][row["met_id2"]] == row["total_score"],
                # check if it is not a self match
                row["met_id1"] != row["met_id2"],
            ]
        ):
            keep.append(i)
    matches = matches.iloc[keep, :]

    if check_charge == True:
        # check whether charge differences have been considered
        keep2 = (
            matches[["charge_diff_inchi", "charge_diff_db", "charge_diff_name"]].sum(
                axis=1
            )
            == 0
        )
        matches = matches.loc[keep2, :]

    # get the pairs and remove the backward matches
    pairs = []
    for i, row in matches.iterrows():
        tpl = (row["met_id1"], row["met_id2"])
        if tpl not in pairs and (tpl[1], tpl[0]) not in pairs:
            pairs.append(tpl)

    return pairs


def removeDuplicateMetabolites(model: MeMoModel) -> tuple[MeMoModel, pd.DataFrame]:
    """Takes a MeMoModel and will remove all the metabolites which are duplicates. Returns the cleaned MeMoModel object"""

    # create a deep copy of the original MeMoMod
    # model = copy.deepcopy(model) # this takes a lot of time and I guess it is not worth it

    # check if the MeMoModel is already annoateted, otherwise do it
    if model.annotated == False:
        model.annotate()
    # write the annotation to the internal cobra model to not loose them
    model.writeAnnotationToCobraModel()

    # do a self match
    matches = model.match(model)

    # get duplicated metabolites
    dups = detectDuplicates(matches)

    # iterate through the metabolites and remove duplicates
    # get all metabolites in the original model - that serves as a check,
    # whether the metabolites have been removed already in a previous iteration
    mets_all_cobra = [x.id for x in model.cobra_model.metabolites]
    change_all = pd.DataFrame()
    for met_id1, met_id2 in dups:
        # get the original ids
        ori_metids1 = [x.orig_ids for x in model.metabolites if x.id == met_id1][0]
        ori_metids2 = [x.orig_ids for x in model.metabolites if x.id == met_id2][0]
        for ori_metid1 in ori_metids1:
            cobra_met1 = model.cobra_model.metabolites.get_by_id(ori_metid1)
            for ori_metid2 in ori_metids2:
                # check if the metabolite is still in the model
                if ori_metid2 in mets_all_cobra:
                    cobra_met2 = model.cobra_model.metabolites.get_by_id(ori_metid2)
                    if cobra_met2.compartment == cobra_met1.compartment:
                        changes_current = mergeAndRemove(cobra_met1, cobra_met2)
                        change_all = pd.concat([change_all, changes_current])
                        # remove the metabolite from the control list
                        mets_all_cobra.remove(ori_metid2)

    memomod = MeMoModel.fromModel(model.cobra_model)

    return (memomod, change_all)
