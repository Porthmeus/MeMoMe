# Porthmeus
# 16.02.24

'''
This file contains function to annotated MeMoMetabolites in a bulk manner. This means a list of MeMoMetabolites is parsed and the annotation takes place on the metabolites in that list. This is mainly to use within the MeMoModel.annotate() function.
'''

import os
import pandas as pd
import warnings
from src.download_db import get_config, get_database_path
from src.MeMoMetabolite import MeMoMetabolite
from src.parseMetaboliteInfos import getAnnoFromIdentifierURL
from src.annotateAux import AnnotationResult


def annotateBiGG_entry(entry:str, database:pd.DataFrame = pd.DataFrame(), allow_missing_dbs: bool = False) -> tuple[dict, list]:
    """
    A small helper function to avoid redundant code.
    Uses a BiGG identifier and annotates it with the identifiers.org entries.
    database - is either empty (default) or a pandas data.frame with the data from the BiGG model homepage for the metabolites. If it is empty, the function will try to load the database from the config file.
    Returns a tuple containing a dictionary for the extracted annotations and a list of new names for the metabolite.
    """

    # check if the database was given, if not, try to load it
    if len(database) == 0:
        # load the database
        config = get_config()
        try:
          db_path =  os.path.join(get_database_path(), config["databases"]["BiGG"]["file"])
          bigg = pd.read_table(db_path)
        except FileNotFoundError as e:
          warnings.warn(str(e))
          # Rethrow exception because we want don't allow missing dbs
          if allow_missing_dbs == False:
            raise e
          return dict(), list()     
    else:
        bigg = database

    annotations = dict()
    # extract the relevant annotation information and return it
    urls = bigg.loc[bigg["universal_bigg_id"]==entry,"database_links"]
    # check if there are entries which are not NA
    urls = urls.loc[~pd.isna(urls)]
    if len(urls) >0:
        urls =";".join(list(urls))
        urls = urls.split(";")
        for url in urls:
            src,val = getAnnoFromIdentifierURL(url)
            if src in annotations.keys():
                annotations[src].append(val)
            else:
                annotations[src] = [val]
    
    # extract the saved name in the database and return it
    names = list(pd.unique(bigg.loc[bigg["universal_bigg_id"]==entry,"name"]))

    return annotations, names
    

def annotateBiGG_id(metabolites: list[MeMoMetabolite], allow_missing_dbs: bool = False) -> AnnotationResult:
    """
    Annotate a list of metabolites with the entries from BiGG. Look for BiGG ids in the metabolite._id slot and if one is found use these. Since BiGG does not provide any InChI strings, this will not results in any new annotated Inchi strings, but it will increase the number of entries in the metabolite annotation.
    """

    # load the database
    config = get_config()
    db_path =  os.path.join(get_database_path(), config["databases"]["BiGG"]["file"])
    try:
        bigg = pd.read_table(db_path)
    except FileNotFoundError as e:
        warnings.warn(str(e))
        # Rethrow exception because we want don't allow missing dbs
        if allow_missing_dbs == False:
          raise e
        return AnnotationResult(0, 0, 0)
    
    new_annos = 0
    new_names = 0
    for met in metabolites:
        if any(bigg["universal_bigg_id"]==met._id):
            new_met_anno_entry,new_names_entry = annotateBiGG_entry(entry = met._id,
                    database = bigg)

            # add new names
            x = met.add_names(new_names_entry)
            new_names = new_names + x

            # add the annotations to the slot in the metabolites
            if len(new_met_anno_entry) > 0:
                x = met.add_annotations(new_met_anno_entry)
                new_annos = new_annos + x 

    # return the number of metabolites which got newly annotated with inchis, annotations and names
    anno_result = AnnotationResult(0, new_annos, new_names)
    return anno_result


def annotateBiGG(metabolites: list[MeMoMetabolite], allow_missing_dbs: bool = False) -> AnnotationResult:
    """
    Annotate a list of metabolites with the entries from BiGG. Look for BiGG ids in the annotation slot of the metabolites and if one is found use these. Since BiGG does not provide any InChI strings, this will not result in any new annotated Inchi strings, but it will increase the number of entries in the metabolite annotation.
    The function will directly add the annotation and names to the MeMoMetabolite object.
    Return value: a tuple of 3 denoting the number of changes metabolites for the following values: [inchi_strings, annotations, names]
    """
    
    # load the database
    config = get_config()
    bigg = None
    try:
      db_path =  os.path.join(get_database_path(), config["databases"]["BiGG"]["file"])
      bigg = pd.read_table(db_path)
    except FileNotFoundError as e:
      warnings.warn(str(e))
      # Rethrow exception because we want don't allow missing dbs
      if allow_missing_dbs == False:
        raise e
      return AnnotationResult(0, 0, 0)
    
    new_annos_added = 0
    new_names_added = 0
    # go through the metabolites and check if there is data which can be added
    for met in metabolites:
        new_met_anno = dict()
        new_names = list() 
        if "bigg.metabolite" in met.annotations.keys():
            for entry in met.annotations["bigg.metabolite"]:
                new_met_anno_entry, new_names_entry = annotateBiGG_entry(entry = entry,
                        database = bigg)
                
                # for all entries create a single dictionary
                for key, value in new_met_anno_entry.items():
                    if key in new_met_anno.keys():
                        new_met_anno[key].extend(value)
                    else:
                        new_met_anno[key] = value
                # combine the names for each entry
                new_names.extend(new_names_entry)

            # add new names to the MeMoMetabolite
            x =met.add_names(new_names)
            new_names_added = new_names_added + x

            # add the annotations to the slot in the metabolites
            if len(new_met_anno) > 0:
                x = met.add_annotations(new_met_anno)
                new_annos_added = new_annos_added + x

    # get the tuple for returning the annotation counts - this is just done to comply to a standard (there will be no InChI annotation from this database
    anno_result = AnnotationResult(0, new_annos_added, new_names_added)
    return anno_result
