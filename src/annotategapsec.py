'''
This file contains function to annotated MeMoMetabolites in a bulk manner. This means a list of MeMoMetabolites is parsed and the annotation takes place on the metabolites in that list. This is mainly to use within the MeMoModel.annotate() function.
'''

import os
import pandas as pd
import warnings
from src.download_db import get_config, get_database_path
from src.MeMoMetabolite import MeMoMetabolite
from src.annotateAux import AnnotationResult
from src.annotateInchiRoutines import findOptimalInchi


def annotategapsec_entry(entry:str, database:pd.DataFrame = pd.DataFrame(), allow_missing_dbs: bool = False) -> list[list, dict, list]:
    """
    A small helper function to avoid redundant code.
    Uses a gapsec identifier and annotates it with the identifiers.org entries.
    database - is either empty (default) or a pandas data.frame with the data from the gapsec homepage for the metabolites. If it is empty, the function will try to load the database from the config file.
    Returns a tuple containing a list of inchi strings, a dictionary for the extracted annotations, and a list of new names for the metabolite.
    """
    # !!! can also extract pka and pkb !!!

    # check if the database was given, if not, try to load it
    if len(database) == 0:
        # load the database
        config = get_config()
        db_path =  os.path.join(get_database_path(), config["databases"]["gapsec"]["file"])
        try:
          gapsec = pd.read_table(db_path)
        except FileNotFoundError as e:
          warnings.warn(e)
          # Rethrow exception because we want don't allow missing dbs
          if allow_missing_dbs == False:
            raise e
          return list(), dict(), list()
    else:
        gapsec = database

    # correct the relevant dataframe column names to conform with the identifiers.org prefixes
    gapsec.rename(columns={"MNX_ID":"metanetx.chemical", "InChIKey":"inchiKey", "hmdbID":"hmdb", "chebiID":"CHEBI", "InChI":"inchi", "keggID":"kegg.compound", "biggID":"bigg.metabolite", "biocycID":"biocyc"}, inplace=True)

    # extract the inchi strings and return them (Note: "InChI" has smaller number of entries "smiles".)
    inchi = gapsec.loc[gapsec["id"]==entry,"inchi"]
    # check if there are entries which are not NA or empty
    inchi = inchi.loc[~pd.isna(inchi)]
    if len(inchi) >0:
       inchi_strings = inchi.to_string(index=False).split(";")
    # get the correct inchi_string, if there was more than one
    if len(inchi_strings) > 1:
      inchi_string = findOptimalInchi(inchi_strings)
      if inchi_string is None:
        raise NotImplementedError()
    elif len(inchi_strings) == 1:
      inchi_string = inchi_strings[0]
    else:
      inchi_string = None

    annotations = dict()
    # extract the relevant annotation information and return it
    keys = ["metanetx.chemical", "inchiKey", "hmdb", "CHEBI", "kegg.compound", "bigg.metabolite", "biocyc"]
    for key in keys:
        values = gapsec.loc[gapsec["id"]==entry,key]
        # check if there are entries which are not NA or empty
        values = values.loc[~pd.isna(values)]
        if len(values) >0:
           values = values.to_string(index=False).split(";")
           if key in annotations.keys():
              annotations[key] += values
           else:
              annotations[key] = values

    # extract the saved name in the database and return it
    names = list(pd.unique(gapsec.loc[gapsec["id"]==entry,"name"]))

    return inchi_string, annotations, names


def annotategapsec_id(metabolites: list[MeMoMetabolite], allow_missing_dbs: bool = False) ->AnnotationResult:
    """
    Annotate a list of metabolites with the entries from gapsec. Look for gapsec ids in the metabolite._id slot and if one is found use these.
    """

    # load the database
    config = get_config()
    db_path =  os.path.join(get_database_path(), config["databases"]["gapsec"]["file"])
    try:
      gapsec = pd.read_table(db_path)
    except FileNotFoundError as e:
      # Rethrow exception because we want don't allow missing dbs
      warnings.warn(str(e))
      if allow_missing_dbs == False:
         raise e
      return AnnotationResult(0, 0, 0)
    
    counter = [0,0,0] # counter for names, annotations, inchi_strings
    for met in metabolites:
        if any(gapsec["id"]==met._id):
            # get the names, annotations, inchi_strings
            new_inchi_string,new_met_annos,new_names = annotategapsec_entry(entry = met._id, 
                                                                            database = gapsec)

            # add the names
            if len(new_names) > 0:
                counter[0] = counter[0] + met.add_names(new_names)

            # add the annotations to the slot in the metabolites
            if len(new_met_annos) > 0:
                counter[1] = counter[1] + met.add_annotations(new_met_annos)

            # add the inchi_string
            if new_inchi_string != None:
                counter[2] = counter[2] + met.add_inchi_string(new_inchi_string)

    anno_result = AnnotationResult(counter[2], counter[1], counter[0])
    return anno_result 


def annotategapsec(metabolites: list[MeMoMetabolite], allow_missing_dbs: bool = False) ->AnnotationResult:
    """
    Annotate a list of metabolites with the entries from gapsec. Look for gapsec ids in the annotation slot of the metabolites and if one is found use these.
    The function will directly add the inchi string, annotations and names to the MeMoMetabolite object.
    Return value: a tuple of 3 denoting the number of changes metabolites for the following values: [inchi string, annotations, names]
     """

    # load the database
    config = get_config()
    db_path =  os.path.join(get_database_path(), config["databases"]["gapsec"]["file"])
    try:
      gapsec = pd.read_table(db_path)
    except FileNotFoundError as e:
      warnings.warn(str(e))
      # Rethrow exception because we want don't allow missing dbs
      if allow_missing_dbs == False:
        raise e
      return AnnotationResult(0, 0, 0)

    counter = [0,0,0] # counter for names, annotations, inchi_strings
    # go through the metabolites and check if there is data which can be added
    for met in metabolites:
        if "seed.compound" in met.annotations.keys():
            met_counter = [0,0,0]
            for entry in met.annotations["seed.compound"]:
                if any(gapsec["id"]==entry):
                   # get the names, annotations, inchi_strings
                   new_inchi_string,new_met_annos,new_names = annotategapsec_entry(entry = entry, 
                                                                             database = gapsec)

                   # add the names
                   if len(new_names) > 0:
                    met_counter[0] = met_counter[0] + met.add_names(new_names)

                   # add the annotations to the slot in the metabolites
                   if len(new_met_annos) > 0:
                    met_counter[1] = met_counter[1] + met.add_annotations(new_met_annos)

                   # add the inchi_string
                   if new_inchi_string != None:
                      met_counter[2] = met_counter[2] + met.add_inchi_string(new_inchi_string)

            # sum the counters
            met_counter = [int(x>0) for x in met_counter]
            counter = [x+y for x,y in zip(counter, met_counter)]

    anno_result = AnnotationResult(counter[2], counter[1], counter[0])
    return anno_result 
