# Porthmeus
# 07.07.23

'''
This file contains function to annotated MeMoMetabolites in a bulk manner. This means a list of MeMoMetabolites is parsed and the annotation takes place on the metabolites in that list. This is mainly to use within the MeMoModel.annotate() function
'''

import numpy as np
import pandas as pd

from src.MeMoMetabolite import MeMoMetabolite
from src.annotateInchiRoutines import findOptimalInchi


def annotateChEBI(metabolites: list[MeMoMetabolite]) -> tuple[int, int]:
    """ Annotate the metaboltes with Inchis from ChEBI """

    # check if any unannotated metabolites exist 
    ids = [x for x, y in enumerate(metabolites) if y._inchi_string == None]
    not_annotated_metabolites = len(ids)
    annotated_by_chebi = 0
    # check if chebis ids are actually present in the annotation slot
    annos = any(["chebi" in x.annotations.keys() for i, x in enumerate(metabolites) if i in ids])
    if annos:
        # download the information from the server
        print("TEST")
        chebi_db = pd.read_table("https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/chebiId_inchi.tsv")
        print("TEST2")
        chebi_db.index = chebi_db['CHEBI_ID']

        # annotate the metabolites with the inchi_string
        for i in ids:
            # For each metabolite that does not have a INCHI string get its chebi id #TODO  KEY OR STRING
            chebis = [int(x.replace("CHEBI:", "")) for x in metabolites[i].annotations["chebi"]]
            # Find all the corresponding INCHI key
            chebis = [x for x in chebis if x in chebi_db.index]
            inchis = np.unique(chebi_db.loc[chebis, "InChI"])
            if len(inchis) > 0:
                inchi = findOptimalInchi(inchis)
                annotated_by_chebi += 1

            metabolites[i].set_inchi_string(inchi)

    return not_annotated_metabolites, annotated_by_chebi


def annotateLove(metabolites: list[MeMoMetabolite]) -> tuple[int, int]:
    """ Annotate the metaboltes with Inchis from ChEBI """

    # check if any unannotated metabolites exist
    ids = [x for x, y in enumerate(metabolites) if y._inchi_string == None]
    not_annotated_metabolites = len(ids)
    print("NOT ANNO", not_annotated_metabolites)

    return 0, 0
