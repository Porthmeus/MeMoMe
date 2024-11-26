'''
This file contains function to annotated MeMoMetabolites in a bulk manner. This means a list of MeMoMetabolites is parsed and the annotation takes place on the metabolites in that list. This is mainly to use within the MeMoModel.annotate() function
'''

import os
import json
import pandas as pd
import warnings
from io import StringIO
from src.download_db import get_config, get_database_path
from src.MeMoMetabolite import MeMoMetabolite
from src.annotation.annotateAux import AnnotationResult, load_database, handleIDs
from typing import Optional

def handle_vmh_entries(vmh, entry):
  keys = ["biggId", "keggId", "cheBlId", "inchiString", "inchiKey", "smile", "hmdb", "metanetx", "seed", "biocyc"]
  vmh = vmh.loc[vmh["abbreviation"] == entry,]
  fullName = list(set(vmh["fullName"]))
  vmh = vmh[keys]
  annotations = dict()
  if len(vmh) > 1:
    raise Exception("More than one value for key: " + entry)

  for index, row in vmh.iterrows():
    for key, value in row.items():
      if (not pd.isna(value)) and (len(str(value)) > 0):
        # If key does not exist, create it and append the value
        annotations.setdefault(key, []).append(value)

  return(annotations, fullName)

def __json_to_tsv(path) -> pd.DataFrame:
  with open(path, "r") as f: 
    vmh = json.load(f)["results"]
    return(pd.DataFrame(vmh))


def annotateVMH_entry(entry: str,  database: pd.DataFrame = pd.DataFrame(), allow_missing_dbs: bool = False) -> tuple[dict, list]:
    """
    A small helper function to avoid redundant code.
    Uses a VMH identifier and annotates it with the identifiers.org entries.
    database - is either empty (default) or a dictionary with the data from the VMH homepage for the metabolites. If it is empty, the function will try to load the database from the config file.
    Returns a tuple containing a dictionary for the extracted annotations and a list of new names for the metabolite.
    """

    # check if the database was given, if not, try to load it
    if len(database) > 0:
      vmh = database
    else:
      vmh = load_database(get_config()["databases"]["VMH"]["file"], allow_missing_dbs, __json_to_tsv)
    
    if vmh.empty:
      return dict(), list()

    annotations, names = handle_vmh_entries(vmh, entry)
    return annotations, names


    



def annotateVMH(metabolites: list[MeMoMetabolite], allow_missing_dbs: bool = False) -> AnnotationResult:
    """
    Annotate a list of metabolites with the entries from VMH. Look for VMH ids in the annotation slot of the metabolites and if one is found use these.
    The function will directly add the annotation and names to the MeMoMetabolite object.
    Return value: a tuple of 3 denoting the number of changes metabolites for the following values: [inchi_strings, annotations, names]
    """
    
    # load the database
    config = get_config()
    db_path =  os.path.join(get_database_path(), config["databases"]["VMH"]["file"])
    try:
      with open(db_path, "r") as f: 
         vmh_json = json.load(f)
    except FileNotFoundError as e:
      warnings.warn(str(e))
      # Rethrow exception because we want don't allow missing dbs
      if allow_missing_dbs == False:
        raise e
      return AnnotationResult(0, 0, 0,)

    vmh = json.dumps(vmh_json['results'])
    vmh = pd.read_json(StringIO(vmh))
    
    new_annos_added = 0
    new_names_added = 0
    # go through the metabolites and check if there is data which can be added
    for met in metabolites:
        new_met_anno = dict()
        new_names = list() 
        if "vmhmetabolite" in met.annotations.keys():
            for entry in met.annotations["vmhmetabolite"]:
                new_met_anno_entry, new_names_entry = annotateVMH_entry(entry = entry,
                        database = vmh_json)
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
    

    # get the tuple for returning the annotation counts
    anno_result = AnnotationResult(0, new_annos_added, new_names_added)
    return anno_result


def annotateVMH_id(metabolites: list[MeMoMetabolite], allow_missing_dbs: bool = False) -> AnnotationResult:
    """
    Annotate a list of metabolites with the entries from VMH. Look for VMH ids in the metabolite._id slot and if one is found use these.
    """

    # check if the database was given, if not, try to load it
    vmh = load_database(get_config()["databases"]["VMH"]["file"], allow_missing_dbs, __json_to_tsv)
    
    if vmh.empty:
      return AnnotationResult(0, 0, 0)
    

    return handleIDs(vmh, metabolites, "abbreviation", annotateVMH_entry)
