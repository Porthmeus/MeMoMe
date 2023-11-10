from src.MeMoModel import *


def identify_remove_rxn_dupl_metid(model, duplcmets, metid2keep=None):

    if metid2keep is None:
        targetid = duplcmets[0]
    else:
        targetid = metid2keep

    rxnDict = {}
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
    forUsrCons = []

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

    return duplcRxn, rxn2remove, forUsrCons


if __name__ == '__main__':
    model = MeMoModel.fromPath(Path(r'C:\Users\Bram Nap\Desktop\Recon3DModel_301.xml'))
    met1 = 'xol7ah3'
    met2 = 'CE0233'
    duplcmets = [met1, met2]
    identify_remove_rxn_dupl_metid(model, duplcmets)
