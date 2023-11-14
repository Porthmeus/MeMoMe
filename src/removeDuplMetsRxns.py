from src.MeMoModel import *
# bramnap 11-2023


def update_model_duplicateMetId(memomemodel, duplcmets, metid2keep=None):
    # TODO: make it so it can take a dictionary with key being the metid2keep?

    """Main script that calls on the modular functions to identify, remove and or update reactions with pre-defined
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
    return upd_rmv_memomemodel


def identify_rxn_dupl_metid(model, duplcmets, targetid):

    metRxn = {}
    easyReadRxnDict = {}
    altEasyReadRxnDict = {}

    for userMet in duplcmets:
        searchID = '^'+userMet+r'\[\w\]'
        modelMets = model.cobra_model.metabolites.query(searchID)

        for met in modelMets:
            rxns = model.cobra_model.metabolites.get_by_id(met.id).reactions
            for rxn in rxns:
                if userMet not in metRxn.keys():
                    metRxn[userMet] = [rxn.id]
                else:
                    metRxn[userMet] += [rxn.id]

                rxninfo = rxn.metabolites
                combinedMets = []
                altcombinedMets = []
                for metstoich in rxninfo:
                    combinedstr = metstoich.id + str(rxninfo[metstoich])
                    altcominedstr = metstoich.id + str(rxninfo[metstoich]*-1)
                    if userMet is not targetid:
                        combinedstr = combinedstr.replace(userMet, targetid)
                        altcominedstr = altcominedstr.replace(userMet, targetid)
                    combinedMets += [combinedstr]
                    altcombinedMets += [altcominedstr]
                    easyReadRxnDict[rxn.id] = sorted(combinedMets)
                    altEasyReadRxnDict[rxn.id] = sorted(altcombinedMets)

    rxns2keep = []
    rxns2inspect = []
    for key in metRxn:
        if key is targetid:
            # Make a dict with all the ids and stoichs per reactions
            rxns2keep += metRxn[key]
        else:
            rxns2inspect += metRxn[key]

    duplcRxn = {}
    rxn2remove = []
    rxn2update = []
    forUsrCons = {}

    for rxnIns in rxns2inspect:
        rxn1 = easyReadRxnDict[rxnIns]
        rxnalt = altEasyReadRxnDict[rxnIns]
        for reactionComp in rxns2keep:
            rxn2 = easyReadRxnDict[reactionComp]
            if rxn1 == rxn2:
                print(rxnIns + ' is the same as ' + reactionComp)

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
                    forUsrCons += [rxnIns]

                if reactionComp not in duplcRxn:
                    duplcRxn[reactionComp] = [rxnIns]
                else:
                    duplcRxn[reactionComp] += [rxnIns]

            elif rxnalt == rxn2:
                print(rxnIns + ' is the same as but reversed ' + reactionComp + ' User has to double check reactions '
                                                                                'before removing')
                forUsrCons += [rxnIns]
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

    for rxn in reactions2remove:
        memome_model.cobra_model.reactions.get_by_id(rxn).upper_bound = 0
        memome_model.cobra_model.reactions.get_by_id(rxn).lower_bound = 0
    return memome_model


def update_reaction(memome_model, reactions_oldMet, newmetID):
    for rxn in reactions_oldMet:
        oldmet = reactions_oldMet[rxn]
                      
        rxn2update = memome_model.cobra_model.reactions.get_by_id(rxn)
        totalmets = rxn2update.metabolites
        oldmetentry = [(key, value) for key, value in totalmets.items() if key.startswith(oldmet)]
        oldmet2remove = oldmetentry[0][0]
        stoich = oldmetentry[0][1]
        comp = oldmet2remove[-3:len(oldmet2remove)]
        met2add = newmetID+comp

        rxn2update.subtract_metabolites({memome_model.cobra_model.metabolites.get_by_id(oldmet2remove): stoich})
        rxn2update.add_metabolites({memome_model.cobra_model.metabolites.get_by_id(met2add): stoich})

    return memome_model


if __name__ == '__main__':
    model = MeMoModel.fromPath(Path(r'C:\Users\Bram Nap\Desktop\Recon3DModel_301.xml'))
    met1 = 'xol7ah3'
    met2 = 'CE0233'
    duplcmets = [met1, met2]
    identify_remove_rxn_dupl_metid(model, duplcmets)
