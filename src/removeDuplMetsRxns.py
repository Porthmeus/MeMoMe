from src.MeMoModel import *
# bramnap 11-2023


def update_model_duplicateMetId(memomemodel, duplcmets, metid2keep=None):
    # TODO: make it so it can take a dictionary with key being the metid2keep?

    """Main function that calls on the modular functions to identify, remove and or update reactions with pre-defined
    duplicated metabolites

    Input:
    memomemodel = MeMoMe model structure. Needs to have the .cobra_model field
    duplcmets = a list of strings with identical metabolites
    metid2keep = string, the metabolite ID to keep in the model

    Output:
    upd_rmv_memomemodel = updated memomemodel
    """

    # If metid2keep is undefined, pick the first element in duplcmets
    if metid2keep is None:
        metid2keep = duplcmets[0]
    else:
        metid2keep = metid2keep

    # Run the different functions to identify, remove and update reactions
    duplcRxn, rxn2remove, forUsrCons, rxn2update = identify_rxn_dupl_metid(memomemodel, duplcmetsm, targetid=metid2keep)
    rmv_memomemodel = remove_reaction(memomemodel, rxn2remove)
    upd_rmv_memomemodel = update_reaction(rmv_memomemodel, xn2update, metid2keep)
    upd_rmv_cobramodel = upd_rmv_memomemodel.cobra_model
    return upd_rmv_memomemodel, upd_rmv_cobramodel


def identify_rxn_dupl_metid(model, duplcmets, targetid):
    """"Function that identifies reactions associated with given metabolites and checks if they're the same reaction
    Input:
    model = MeMoMe model structure. Needs to have the .cobra_model field
    duplcmets = a list of strings with identical metabolites
    targetid = string, the metabolite ID to keep in the model

    Output:
    duplcRxn = dict, stores which reactions are the same {rxn desired ID: [rxn(s) undesired ID]}
    rxn2remove = list, reactions that have to be removed as they are duplicates
    forUsrCons = dict, duplicated reactions that need to be checked for bound or gpr information before removing
                 {rxn desired ID: [rxn(s) undesired ID]}
    rxn2update = dict, unique reactions that need to have the metabolite ID updated {rxn ID: undesired met ID}
    """

    # Initialise storage variables
    metRxn = {}
    easyReadRxnDict = {}
    altEasyReadRxnDict = {}

    # Find and save in which compartments the metabolites are present
    for userMet in duplcmets:
        searchID = '^'+userMet+r'\[\w\]'
        modelMets = model.cobra_model.metabolites.query(searchID)

        # Obtain the reactions associated with compartment specific metabolites
        for met in modelMets:
            rxns = model.cobra_model.metabolites.get_by_id(met.id).reactions
            for rxn in rxns:
                if userMet not in metRxn.keys():
                    metRxn[userMet] = [rxn.id]
                else:
                    metRxn[userMet] += [rxn.id]

                # For associated reactions obtain all metabolites and stoichiometries that play a role
                rxninfo = rxn.metabolites
                combinedMets = []
                altcombinedMets = []
                for metstoich in rxninfo:
                    combinedstr = metstoich.id + str(rxninfo[metstoich])

                    # Account for eventual reversal of substrates and products
                    altcominedstr = metstoich.id + str(rxninfo[metstoich]*-1)
                    if userMet is not targetid:
                        combinedstr = combinedstr.replace(userMet, targetid)
                        altcominedstr = altcominedstr.replace(userMet, targetid)
                    combinedMets += [combinedstr]
                    altcombinedMets += [altcominedstr]
                    easyReadRxnDict[rxn.id] = sorted(combinedMets)
                    altEasyReadRxnDict[rxn.id] = sorted(altcombinedMets)

    # Initialse variables
    rxns2keep = []
    rxns2inspect = []

    # Seperate reactions that have the desired and undesired metabolite IDs
    for key in metRxn:
        if key is targetid:
            # Make a dict with all the ids and stoichs per reactions
            rxns2keep += metRxn[key]
        else:
            rxns2inspect += metRxn[key]

    # Initialise variables
    duplcRxn = {}
    rxn2remove = []
    rxn2update = {}
    forUsrCons = {}

    # Check if the two reactions are the same metabolite and stoichiometry wise
    for rxnIns in rxns2inspect:
        rxn1 = easyReadRxnDict[rxnIns]
        rxnalt = altEasyReadRxnDict[rxnIns]
        for reactionComp in rxns2keep:
            rxn2 = easyReadRxnDict[reactionComp]
            if rxn1 == rxn2:
                print(rxnIns + ' is the same as ' + reactionComp)

                # Check additional information on reactions to ensure no information is lost
                rxn1ub = model.cobra_model.reactions.get_by_id(rxnIns).upper_bound
                rxn1lb = model.cobra_model.reactions.get_by_id(rxnIns).lower_bound
                rxn1gpr = model.cobra_model.reactions.get_by_id(rxnIns).gene_reaction_rule

                rxn2ub = model.cobra_model.reactions.get_by_id(reactionComp).upper_bound
                rxn2lb = model.cobra_model.reactions.get_by_id(reactionComp).lower_bound
                rxn2gpr = model.cobra_model.reactions.get_by_id(reactionComp).gene_reaction_rule

                if rxn1ub == rxn2ub and rxn1lb == rxn2lb and rxn1gpr == rxn2gpr:
                    rxn2remove += [rxnIns]
                else:
                    print(f'The bounds and/or gpr for {rxnIns} and {reactionComp} are not the same. User should check the '
                          f'reactions and make sure all information is stored in 1 reaction before removing the other')
                    if reactionComp not in forUsrCons:
                        forUsrCons[reactionComp] = [rxnIns]
                    else:
                        forUsrCons[reactionComp] += [rxnIns]

                if reactionComp not in duplcRxn:
                    duplcRxn[reactionComp] = [rxnIns]
                else:
                    duplcRxn[reactionComp] += [rxnIns]

            elif rxnalt == rxn2:
                print(rxnIns + ' is the same as but reversed ' + reactionComp + ' User has to double check reactions '
                                                                                'before removing')
                if reactionComp not in forUsrCons:
                    forUsrCons[reactionComp] = [rxnIns]
                else:
                    forUsrCons[reactionComp] += [rxnIns]

                if reactionComp not in duplcRxn:
                    duplcRxn[reactionComp] = [rxnIns]
                else:
                    duplcRxn[reactionComp] += [rxnIns]
            else:
                for key in metRxn:
                    if rxnIns in metRxn[key]:
                        rxn2update[rxnIns] = [key]

    return duplcRxn, rxn2remove, forUsrCons, rxn2update


def remove_reaction(memome_model, reactions2remove):
    """" Function that "removes" reactions - now only sets the bounds to zero
    Input:
    memome_model = MeMoMe model structure. Needs to have the .cobra_model field
    reactions2remove = list, reaction IDs to be removed

    Output:
    memome_model = MeMoMe model structure. updated to have reactions2remove bounds set to 0
    """

    for rxn in reactions2remove:
        memome_model.cobra_model.reactions.get_by_id(rxn).upper_bound = 0
        memome_model.cobra_model.reactions.get_by_id(rxn).lower_bound = 0
    return memome_model


def update_reaction(memome_model, reactions_oldMet, newmetID):
    """" Function that updated the metabolite ID in specified reactions
    Input:
    memome_model = MeMoMe model structure. Needs to have the .cobra_model field
    reactions_oldMet = dict, reactions that need to have the metabolite ID updated {rxn ID: undesired met ID}
    newmetID = string, desired metabolite ID

    Output:
    memome_model = MeMoMe model structure. updated have specified reactions carry the desired metabolite ID
    """

    for rxn in reactions_oldMet:
        # Obtain old metabolite ID
        oldmet = reactions_oldMet[rxn]

        # set the reaction object from the model for easy use
        rxn2update = memome_model.cobra_model.reactions.get_by_id(rxn)
        # Obtain the metabolites of the reaction and the comparment and stoichiometry of the undesired metabolite ID
        totalmets = rxn2update.metabolites
        oldmetentry = [(key, value) for key, value in totalmets.items() if key.startswith(oldmet)]
        oldmet2remove = oldmetentry[0][0]
        stoich = oldmetentry[0][1]
        comp = oldmet2remove[-3:len(oldmet2remove)]
        # Set the compartment to the desired metabolite ID
        met2add = newmetID+comp

        # Remove and add the metabolites to the reactions
        rxn2update.subtract_metabolites({memome_model.cobra_model.metabolites.get_by_id(oldmet2remove): stoich})
        rxn2update.add_metabolites({memome_model.cobra_model.metabolites.get_by_id(met2add): stoich})

    return memome_model


if __name__ == '__main__':
    model = MeMoModel.fromPath(Path(r'C:\Users\Bram Nap\Desktop\Recon3DModel_301.xml'))
    met1 = 'xol7ah3'
    met2 = 'CE0233'
    duplcmets = [met1, met2]
    identify_remove_rxn_dupl_metid(model, duplcmets)
