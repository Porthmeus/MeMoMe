# Porthmeus
# 07.07.23

'''
This file contains function to annotated MeMoMetabolites in a bulk manner. This means a list of MeMoMetabolites is parsed and the annotation takes place on the metabolites in that list. This is mainly to use within the MeMoModel.annotate() function
'''

import numpy as np
import pandas as pd
import logging
from src.MeMoMetabolite import MeMoMetabolite
from src.annotateInchiRoutines import findOptimalInchi

logger = logging.getLogger('logger')


def annotateChEBI(metabolites: list[MeMoMetabolite]) -> tuple[int, int]:
    """ Annotate the metaboltes with Inchis from ChEBI """

    # check if any unannotated metabolites exist 
    ids = [x for x, y in enumerate(metabolites) if y._inchi_string == None]
    not_annotated_metabolites = len(ids)
    annotated_by_chebi = 0
    # check if chebis ids are actually present in the annotation slot
    # TODO This does not work with 35, because we try to access the "chebi" even for thos metabolies that don;t have one
    annos = any(["chebi" in x.annotations.keys() for i, x in enumerate(metabolites) if i in ids])
    if annos:
        # download the information from the server
        chebi_db = pd.read_table("/home/td/Downloads/chebiId_inchi.tsv")
        chebi_db.index = chebi_db['CHEBI_ID']

        # annotate the metabolites with the inchi_string
        # For each metabolite that does not have a INCHI string get its chebi id
        for i in ids:
            # If the current metabolite does not have an chebi annotation skip it
            if "chebi" not in metabolites[i].annotations:
                continue
            chebis = [int(x.replace("CHEBI:", "")) for x in metabolites[i].annotations["chebi"]]
            # Find all the corresponding INCHI key
            chebis = [x for x in chebis if x in chebi_db.index]
            inchis = np.unique(chebi_db.loc[chebis, "InChI"])
            if len(inchis) > 0:
                inchi = findOptimalInchi(inchis)
                annotated_by_chebi += 1

            metabolites[i].set_inchi_string(inchi)

    return not_annotated_metabolites, annotated_by_chebi


def annotateVMH_HMDB(metabolites: list[MeMoMetabolite]) -> tuple[int, int]:
    """ Annotate the metaboltes with Inchis from """
    vmh_db = pd.read_json("/home/td/Projects/MeMoMe/Databases/vmh.json")
    res = vmh_db["results"]
    vmh_db = pd.DataFrame.from_records(res)

    # check if any unannotated metabolites exist
    ids = [x for x, y in enumerate(metabolites) if y._inchi_string == None]
    not_annotated_metabolites = len(ids)
    print("Not annotated before VMH " + str(not_annotated_metabolites))
    annotate_counter_hmdb = 0
    annos = any(["hmdb" in x.annotations.keys() for i, x in enumerate(metabolites) if i in ids])

    for i in ids:
        #####################
        #
        #                       HMDB
        #
        #####################
        # If the current metabolite does not have a hmdb annotation skip it
        if "hmdb" not in metabolites[i].annotations:
            continue
        # Get the bigg ids for each metabolite
        seeds = metabolites[i].annotations["hmdb"]

        # TODO CHECK FOR UNIQUE LENGTH
        if len(seeds) > 1:
            print("THIS DOES HAPPEN: WUHU")
        id_ = seeds[0]
        matches: pd.DataFrame = vmh_db[vmh_db['hmdb'] == id_]
        if len(matches.index) > 0:
            # Cool we got a match, so we don't have to try the longer id
            # But we still have to figure out if there is an INCHI string associated to it
            inchi_strings = matches["inchiString"]
            inchi_strings = list(inchi_strings)
            inchi_strings.remove('')
            if len(inchi_strings) > 0:
                metabolites[i].set_inchi_string(findOptimalInchi(inchi_strings))
                annotate_counter_hmdb += 1
                continue
        else:
            splits = id_.split("B")
            assert len(splits) == 2
            missing_zeros = 7 - len(splits[1])
            new_id = splits[0] + "B" + "0" * missing_zeros + splits[1]
            matches: pd.DataFrame = vmh_db[vmh_db['hmdb'] == new_id]
            if len(matches.index) > 0:
                # Cool we got a match with the new id so we check if it contains an INCHI string
                # If yes, we continue with the next metabolite otherwise.
                # If no, we have to try the next useful identifier
                inchi_strings = matches["inchiString"]
                inchi_strings = list(inchi_strings)
                if '' in inchi_strings:
                    inchi_strings.remove('')
                if len(inchi_strings) > 0:
                    metabolites[i].set_inchi_string(findOptimalInchi(inchi_strings))
                    annotate_counter_hmdb += 1
                    continue
                pass
    logger.info("Annotated by HMDB " + str(annotate_counter_hmdb))

    ids = [x for x, y in enumerate(metabolites) if y._inchi_string == None]
    return not_annotated_metabolites, annotate_counter_hmdb


