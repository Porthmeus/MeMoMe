# Porthmeus
# 07.07.23

'''
This file contains function to annotated MeMoMetabolites in a bulk manner. This means a list of MeMoMetabolites is parsed and the annotation takes place on the metabolites in that list. This is mainly to use within the MeMoModel.annotate() function
'''

import os
import numpy as np
import pandas as pd
from pathlib import Path
from src.download_db import get_config, get_database_path



from src.MeMoMetabolite import MeMoMetabolite
from src.annotateInchiRoutines import findOptimalInchi
from src.download_db import databases_available, get_config, get_database_path
from src.parseMetaboliteInfos import getAnnoFromIdentifierURL



def annotateChEBI(metabolites: list[MeMoMetabolite]) -> tuple[int, int]:
    """ Annotate the metaboltes with Inchis from ChEBI """

    # check if databases are available and if so load the config
    if databases_available():
        config = get_config()
        db_path = get_database_path()

    # check if any unannotated metabolites exist and whether these have a ChEBI entry
    ids = [x for x, y in enumerate(metabolites) if y._inchi_string == None and "chebi" in y.annotations.keys()]
    not_annotated_metabolites = len(ids)
    annotated_by_chebi = 0
    # check if chebis ids are actually present in the annotation slot
    # TODO This does not work with 35, because we try to access the "chebi" even for thos metabolies that don;t have one
    annos = any(["chebi" in x.annotations.keys() for i, x in enumerate(metabolites) if i in ids])
    if annos:
        # download the information from the server
        #chebi_db = pd.read_table("/home/td/Downloads/chebiId_inchi.tsv")
        chebi_db = pd.read_table(db_path.joinpath(config["databases"]["ChEBI"]["file"]))
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

def annotateBiGG_entry(entry:str, annotations:dict = dict(), database:pd.DataFrame = pd.DataFrame()) -> dict:
    """
    A small helper function to avoid redundant code
    Uses a BiGG identifiers and annotates it with the identifier.org entries.
    The results are stored int he annotations dictionary. This is empty by
    default, but if it exists, the entries will be simply added (duplicates are
    allowed here).
    database - is either empty (default) or a pandas data.frame with the data from the BiGG model homepage for the metabolites. If it is empty, the function will try to load the database from the config file
    """
    # check if the database was given, if not, try to load it
    if len(database) == 0:
        # load the database
        config = get_config()
        db_path =  os.path.join(get_database_path(), config["databases"]["BiGG"]["file"])
        bigg = pd.read_table(db_path)
    else:
        bigg = database

    # extract the relevant annotation information and return it
    urls = bigg.loc[bigg["universal_bigg_id"]==entry,"database_links"]
    # check if there are entries which are not NA
    urls = urls.loc[~pd.isna(urls)]
    if len(urls) >0:
        urls =";".join(list(urls))
        urls = urls.split(";")
        for url in urls:
            src,val =getAnnoFromIdentifierURL(url)
            if src in annotations.keys():
                annotations[src].append(val)
            else:
                annotations[src] = [val]
    
    # extract the saved name in the database and return it
    names = list(pd.unique(bigg.loc[bigg["universal_bigg_id"]==entry,"name"]))
    return(annotations, names)

    
    

def annotateBiGG(metabolites: list[MeMoMetabolite]) -> tuple[int, int]:
    """
    Annotate a list of metabolites with the entries from BiGG. Look for BiGG ids in the annotation slot of the metabolites and if one is found use these. Since BiGG does not provide any InChI strings, this will not results in any new annotated Inchi strings, but it will increase the number of entries in the metabolite annotation.
    """
    
    # load the database
    config = get_config()
    db_path =  os.path.join(get_database_path(), config["databases"]["BiGG"]["file"])
    bigg = pd.read_table(db_path)
    
    counter = 0
    # go through the metabolites and check if there is data which can be added
    for met in metabolites:
        new_met_anno = dict()
        if "bigg.metabolite" in met.annotations.keys():
            for entry in met.annotations["bigg.metabolite"]:
                new_met_anno, new_names = annotateBiGG_entry(entry = entry,
                        annotations = new_met_anno,
                        database = bigg)
                # add new names
                met.add_names(new_names)

            # remove duplicates from the annotations
            for key in new_met_anno.keys():
               new_val = list(set(new_met_anno[key]))
               new_met_anno[key] = new_val
            
            # add the annotations to the slot in the metabolites
            if len(new_met_anno) > 0:
                met.add_annotations(new_met_anno)
                counter = counter +1
    
    # return a message how many annotations could be found
    print("Found " + str(counter) + " annotations in the BiGG database, which have been added")

    # get the tuple for returning the annotation counts - this is just done to comply to a standard (there will be no InChI annotation from this database
    ret1,ret2 = annotateLove(metabolites)
    return ret1,ret2


def annotateBiGG_id(metabolites: list[MeMoMetabolite]) -> tuple[int, int]:
    """
    Annotate a list of metabolites with the entries from BiGG. Look for BiGG ids in the metabolite._id slot and if one is found use these. Since BiGG does not provide any InChI strings, this will not results in any new annotated Inchi strings, but it will increase the number of entries in the metabolite annotation.
    """

    # load the database
    config = get_config()
    db_path =  os.path.join(get_database_path(), config["databases"]["BiGG"]["file"])
    bigg = pd.read_table(db_path)
    
    counter = 0
    for met in metabolites:
        new_met_anno = dict()
        #if met._id in bigg["universal_bigg_id"]: # does not work, for whatever reason
        if any(bigg["universal_bigg_id"]==met._id):
            new_met_anno,new_names = annotateBiGG_entry(entry = met._id,
                    annotations = new_met_anno,
                    database = bigg)
            # add new names
            met.add_names(new_names)

            # remove duplicates from the annotations
            for key in new_met_anno.keys():
               new_val = list(set(new_met_anno[key]))
               new_met_anno[key] = new_val
        
            # add the annotations to the slot in the metabolites
            if len(new_met_anno) > 0:
                met.add_annotations(new_met_anno)
                counter = counter +1

    # return a message how many annotations could be found
    print("Found " + str(counter) + " annotations in the BiGG database, which have been added")

    # get the tuple for returning the annotation counts - this is just done to comply to a standard (there will be no InChI annotation from this database
    ret1,ret2 = annotateLove(metabolites)
    return ret1,ret2
