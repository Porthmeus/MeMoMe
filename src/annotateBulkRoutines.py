# Porthmeus
# 07.07.23

'''
This file contains function to annotated MeMoMetabolites in a bulk manner. This means a list of MeMoMetabolites is parsed and the annotation takes place on the metabolites in that list. This is mainly to use within the MeMoModel.annotate() function
'''
import logging

import numpy as np
import pandas as pd

from src.MeMoMetabolite import MeMoMetabolite
from src.annotateInchiRoutines import findOptimalInchi

logger = logging.getLogger('logger')


def annotateChEBI(metabolites: list[MeMoMetabolite]) -> tuple[int, int]:
    """ Annotate the metaboltes with Inchis from ChEBI """

    # check if any unannotated metabolites exist 
    ids = [x for x, y in enumerate(metabolites) if y._inchi_string == None]
    if len(ids) == 0:
        return len(ids), 0
    not_annotated_metabolites = len(ids)
    logger.info("Not annotated before PBC " + str(not_annotated_metabolites))

    annotated_by_chebi = 0
    # check if chebis ids are actually present in the annotation slot
    # TODO This does not work with 35, because we try to access the "chebi" even for thos metabolies that don;t have one
    annos = any(["chebi" in x.annotations.keys() for i, x in enumerate(metabolites) if i in ids])
    if annos:
        # download the information from the server
        #chebi_db = pd.read_table("/home/td/Downloads/chebiId_inchi.tsv")
        chebi_db = pd.read_table("https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/chebiId_inchi.tsv")
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
    vmh_db = pd.read_json("../Databases/vmh.json")
    res = vmh_db["results"]
    vmh_db = pd.DataFrame.from_records(res)

    # check if any unannotated metabolites exist
    ids = [x for x, y in enumerate(metabolites) if y._inchi_string == None]
    if len(ids) == 0:
        return len(ids), 0
    not_annotated_metabolites = len(ids)
    logger.info("Not annotated before VMH " + str(not_annotated_metabolites))
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
    return not_annotated_metabolites, annotate_counter_hmdb


def annotate_PBC(metabolites: list[MeMoMetabolite]) -> tuple[int, int]:
    """ Annotate the metaboltes with Inchis from """
    # check if any unannotated metabolites exist
    # Apparently this list contains duplicates 9/10 are both 10-hydrofolate sth in Dafnis Models
    ids = [x for x, y in enumerate(metabolites) if y._inchi_string == None]
    if len(ids) == 0:
        return len(ids), 0
    # Only load DB if we really need it
    pbc_db = pd.read_csv("Databases/all_compounds.csv", header = None)
    not_annotated_metabolites = len(ids)
    logger.info("Not annotated before PBC " + str(not_annotated_metabolites))
    annotate_counter_pbc = 0
    have_pub_but_no_inchi = ["pubchem.compound" in x.annotations.keys() for i, x in enumerate(metabolites) if i in ids]
    # Just for testing purpose, tells use how many metabolites do not have an inchi annotation but do have a pubmed id
    c = have_pub_but_no_inchi.count(True)
    logger.info("Amount of keys that have a PBC id but no inchi " + str(c))
    annos = any(have_pub_but_no_inchi)

    for i in ids:
        #####################
        #
        #                       Pubchem
        #
        #####################
        # If the current metabolite does not have a PUBCHEM annotation skip it
        if "pubchem.compound" not in metabolites[i].annotations:
            continue
        # Get the PUBHCEM ids for each metabolite
        seeds = metabolites[i].annotations["pubchem.compound"]
        if len(seeds) > 1:
            print("This metabolite has more then one PBC id. Please handle. FOR NOW I will ignore this and take the first id")
        id_ = int(seeds[0])
        matches: pd.DataFrame = pbc_db[pbc_db.iloc[:, 0] == id_]
        if len(matches) > 1:
            raise NotImplementedError("UPSIDUS")

        inchi_strings = matches.iloc[:, 1]

        if len(inchi_strings) > 0:
            metabolites[i].set_inchi_string(findOptimalInchi(inchi_strings))
            annotate_counter_pbc += 1
            continue
        else:
            # Metabolites that have a PBC id but no inchi
            print(seeds)

    return not_annotated_metabolites,annotate_counter_pbc


