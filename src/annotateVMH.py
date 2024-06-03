'''
This file contains function to annotated MeMoMetabolites in a bulk manner. This means a list of MeMoMetabolites is parsed and the annotation takes place on the metabolites in that list. This is mainly to use within the MeMoModel.annotate() function.
'''

import os
import json
import pandas as pd
import warnings
from io import StringIO
from src.download_db import get_config, get_database_path
from src.MeMoMetabolite import MeMoMetabolite
from src.annotateAux import AnnotationResult
from src.annotateInchiRoutines import findOptimalInchi


def correctAnnotationKeys(data:dict) -> dict:
    """
    A small helper function.
    Takes the VMH data dictionary and corrects the relevant keys in the dictionary to conform with the identifiers.org prefixes.
    """

    key_map = {"biggId":"bigg.metabolite", "keggId":"kegg.compound", "cheBlId":"CHEBI", "inchiString":"inchi", "metanetx":"metanetx.chemical", "seed":"seed.compound"}
    new_data = dict()
    for key, value in data.items():
       new_key = key_map.get(key, key)
       new_data[new_key] = value

    return(new_data)


def annotateVMH_entry(entry:str, database:dict = dict(), allow_missing_dbs: bool = False) -> tuple[list, dict, list]:
    """
    A small helper function to avoid redundant code.
    Uses a VMH identifier and annotates it with the identifiers.org entries.
    database - is either empty (default) or a dictionary with the data from the VMH homepage for the metabolites. If it is empty, the function will try to load the database from the config file.
    Returns a tuple containing a list of inchi strings, a dictionary for the extracted annotations, and a list of names for the metabolite.
    """

    # check if the database was given, if not, try to load it
    if len(database) == 0:
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
          return list(), dict(), list()
        vmh = vmh_json['results']
    else:
        vmh_json = database
        vmh = vmh_json['results']

    # find the metabolite information in the database
    for index, dictionary in enumerate(vmh):
        if dictionary["abbreviation"]==entry:
            data = vmh[index]
            new_data = correctAnnotationKeys(data)

    inchi_strings = []
    # extract the inchi strings and return them (Note: "inchiString" and "smile" have the same number of entries.)
    inchi = new_data["inchi"]
    # check if there are entries which are not NA or empty
    if (not pd.isna(inchi)) and (len(str(inchi)) >0):
       inchi_strings.append(inchi)
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
    keys = ["bigg.metabolite", "kegg.compound", "CHEBI", "inchiKey", "hmdb", "metanetx.chemical", "seed.compound", "biocyc"]
    for key in keys:
        value = new_data[key]
        # check if there are entries which are not NA or empty
        if (not pd.isna(value)) and (len(str(value)) >0):
            if key in annotations.keys():
                annotations[key].append(value)
            else:
                annotations[key] = [value]

    # extract the saved name in the database and return it
    names = list(set(new_data["fullName"].split()))

    return inchi_string, annotations, names


def annotateVMH_id(metabolites: list[MeMoMetabolite], allow_missing_dbs: bool = False) -> AnnotationResult:
    """
    Annotate a list of metabolites with the entries from VMH. Look for VMH ids in the metabolite._id slot and if one is found use these.
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
      return AnnotationResult(0, 0, 0)

    vmh = json.dumps(vmh_json['results'])
    vmh = pd.read_json(StringIO(vmh))
    
    counter = [0,0,0] # counter for names, annotations, inchi_strings
    for met in metabolites:
        if any(vmh["abbreviation"]==met._id):
            # get the names, annotations, inchi_strings
            new_inchi_string,new_met_annos,new_names = annotateVMH_entry(entry = met._id,
                    database = vmh_json)

            # add the names
            if len(new_names) > 0:
                counter[0] = counter[0] + met.add_names(new_names)

            # add the annotations to the slot in the metabolites
            if len(new_met_annos) > 0:
                counter[1] = counter[1] + met.add_annotations(new_met_annos)

            # add the inchi_string
            if new_inchi_string != None:
                counter[2] = counter[2] + met.add_inchi_string(new_inchi_string)

    # return the number of metabolites which got newly annotated with inchis, annotations and names
    anno_result = AnnotationResult(counter[2], counter[1], counter[0])
    return anno_result


def annotateVMH(metabolites: list[MeMoMetabolite], allow_missing_dbs: bool = False) -> AnnotationResult:
    """
    Annotate a list of metabolites with the entries from VMH. Look for VMH ids in the annotation slot of the metabolites and if one is found use these.
    The function will directly add the annotation and names to the MeMoMetabolite object.
    Return value: a tuple of 3 denoting the number of changes metabolites for the following values: [inchi string, annotations, names]
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
      return AnnotationResult(0, 0, 0)

    vmh = json.dumps(vmh_json['results'])
    vmh = pd.read_json(StringIO(vmh))
    
    counter = [0,0,0] # counter for names, annotations, inchi_strings
    # go through the metabolites and check if there is data which can be added
    for met in metabolites:
        if "vmhmetabolite" in met.annotations.keys():
            met_counter = [0,0,0]
            for entry in met.annotations["vmhmetabolite"]:
                if any(vmh["abbreviation"]==entry):
                   # get the names, annotations, inchi_strings
                   new_inchi_string,new_met_annos,new_names = annotateVMH_entry(entry = entry, 
                                                                             database = vmh_json)

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

    # get the tuple for returning the annotation counts
    anno_result = AnnotationResult(counter[2], counter[1], counter[0])
    return anno_result
