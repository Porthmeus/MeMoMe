# Porthmeus
# 28.05.24

import re
import cobra as cb
import copy

def printModelStats(model:cb.Model):
    print(f'{len(model.reactions)} reactions')
    print(f'{len(model.metabolites)} metabolites')
    print(f'{len(model.genes)} genes')

def findCommonReactions(met1:cb.Metabolite, met2:cb.Metabolite, reversible_is_same = True) -> list[cb.Reaction]:
    '''Takes two metabolites of the same model and returns reactions which are the same except for the two metabolites
    met1,met2 - the two cobra metabolite objects which are the categorized as same
    reverse_is_same - boolean, should reactions be reported if one of them is reversible but the other is not'''

    # check if the metabolites originate from the same model
    if met1.model != met2.model:
        raise ValueError("Metabolites do not origin from the same model! Aborted.")
    
    # get the essential data
    reactions1 = [x for x in met1.reactions]
    reactions2 = [x for x in met2.reactions]
    reaction_rules1 = [x.reaction for x in reactions1]
    reaction_rules2 = [x.reaction for x in reactions2]
    metId1 = met1.id
    metId2 = met2.id
    
    # iterate through list of reactions
    same =[]
    for i in range(len(reactions1)):
        for j in range(len(reactions2)):
            reac_rule1 = copy.deepcopy(reaction_rules1[i])
            reac_rule2 = copy.deepcopy(reaction_rules2[j])
            if reversible_is_same==True:
                if reac_rule1.find("<=>") != -1 or reac_rule2.find("<=>") != -1:
                    reac_rule1 = reac_rule1.replace("-->","<=>").replace("<--","<=>")
                    reac_rule2 = reac_rule2.replace("-->","<=>").replace("<--","<=>")
            else:
                print(reac_rule1)
                print(reac_rule2)
            if reac_rule1 == re.sub(metId2, metId1, reac_rule2):
                same.append((reactions1[i], reactions2[j]))
    return(same)

def mergeReactions(rxn1:cb.Reaction, rxn2:cb.Reaction, remove_orphans:bool = True, verbose = False) -> None:
    ''' Takes a pair of reactions from the same model and merges them. This will result in merging the information into the first reaction and deletion of the second.
    
    remove_orphans = boolean, if True (default) - all orphanized metabolites and genes in the model will be deleted
    '''
    
    # check if compartments are the same
    if len(rxn1.compartments) > 1:
        raise NotImplementedError("Reaction "+ rxn1.id + " is associated to more than one compartment, don't know how to handle that")
    if len(rxn2.compartments) > 1:
        raise NotImplementedError("Reaction "+ rxn2.id + " is associated to more than one compartment, don't know how to handle that")
    assert (rxn1.compartments == rxn2.compartments), f"You try to merge reactions which take place in different compartments ({rxn1.id=}, {rxn2.id =})"

    
    # check if the gprs are the same, if not merge
    if rxn1.gene_reaction_rule != rxn2.gene_reaction_rule:
        rxn1.gene_reaction_rule = "(" + rxn1.gene_reaction_rule + ") or (" + rxn2.gene_reaction_rule + ")"

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




def exchangeMetabolite(met1:cb.Metabolite, met2:cb.Metabolite) -> None:
    '''Find the all reactions which use the second metabolite and exchange it with the first metabolite'''
    rxns = met2.reactions
    for rxn in rxns:
        coef = rxn.metabolites[met2]
        sub_dic = {met2 : coef,
                met1 : -1*coef}
        rxn.subtract_metabolites(sub_dic)

