# Porthmeus
# 16.02.24

'''
This file contains function to annotated MeMoMetabolites in a bulk manner. This means a list of MeMoMetabolites is parsed and the annotation takes place on the metabolites in that list. This is mainly to use within the MeMoModel.annotate() function
'''

import numpy as np
import pandas as pd
import warnings
from src.MeMoMetabolite import MeMoMetabolite
from src.annotation.annotateInchiRoutines import findOptimalInchi
from src.download_db import get_config, get_database_path
from src.annotation.annotateAux import AnnotationResult, load_database

def annotateChEBI(metabolites: list[MeMoMetabolite], allow_missing_dbs: bool = False) -> AnnotationResult:
    """ Annotate the metaboltes with Inchis from ChEBI """

    chebi_db =  load_database(get_config()["databases"]["ChEBI"]["file"], 
                          allow_missing_dbs, 
                          lambda path: pd.read_csv(path, sep="\t"))
    if chebi_db.empty:
      return AnnotationResult(0, 0,0 )

    # check if any unannotated metabolites exist and whether these have a ChEBI entry
    ids = [x for x, y in enumerate(metabolites) if y._inchi_string == None and "chebi" in y.annotations.keys()]
    not_annotated_metabolites = len(ids)
    annotated_by_chebi = 0
    # check if chebis ids are actually present in the annotation slot
    annos = any(["chebi" in x.annotations.keys() for i, x in enumerate(metabolites) if i in ids])
    if annos:
        chebi_db.index = chebi_db['CHEBI_ID']

        # annotate the metabolites with the inchi_string
        for i in ids:
            # For each metabolite that does not have a INCHI string get its chebi id #TODO  KEY OR STRING
            chebis = [int(x.replace("CHEBI:", "")) for x in metabolites[i].annotations["chebi"]]
            # Find all the corresponding INCHI key
            chebis = [x for x in chebis if x in chebi_db.index]
            inchis = np.unique(chebi_db.loc[chebis, "InChI"])
            if len(inchis) > 0:
                inchi = findOptimalInchi(inchis, charge = metabolites[i]._charge)
                if inchi is None:
                    continue
            else :
                inchi = None

            annotated_by_chebi = annotated_by_chebi + metabolites[i].set_inchi_string(inchi, source = "chebi")
    
    anno_result = AnnotationResult(annotated_by_chebi, 0, 0)
    return anno_result 
