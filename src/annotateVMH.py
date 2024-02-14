import os
import json
import pandas as pd
from io import StringIO
from src.download_db import get_config, get_database_path
from src.MeMoMetabolite import MeMoMetabolite


def annotateVMH_entry(entry:str, annotations:dict = dict(), database:pd.DataFrame = pd.DataFrame()) -> dict:
    """
    A small helper function to avoid redundant code
    Uses a VMH identifier and annotates it with the identifiers.org entries.
    The results are stored in the annotations dictionary. This is empty by
    default, but if it exists, the entries will be simply added (duplicates are
    allowed here).
    database - It is either empty (default) or a pandas data.frame with the data from the VMH homepage for the metabolites. If it is empty, the function will try to load the database from the config file.
    """
    
    # check if the database was given, if not, try to load it
    if len(database) == 0:
        # load the database
        config = get_config()
        db_path =  os.path.join(get_database_path(), config["databases"]["VMH"]["file"])
        with open(db_path) as f:
            vmh_json = json.load(f)
        vmh_json = json.dumps(vmh_json['results'])
        vmh = pd.read_json(StringIO(vmh_json))
    else:
        vmh = database

    # extract the relevant annotation information and return it
    keys = ["biggId", "keggId", "cheBlId", "inchiString", "inchiKey", "smile", "hmdb", "metanetx", "seed", "biocyc"]
    for key in keys:
        k = vmh.loc[vmh["abbreviation"]==entry,key]
        # check if there are entries which are not NA
        k = k.loc[~pd.isna(k)]
        if len(k) >0:
            value = k.values[0] # or value = ("{0}".format("', '".join(k.values)))
            if len(value) >0:
                if key in annotations.keys():
                    annotations[key].append(value)
                else:
                    annotations[key] = [value]

    # extract the saved name in the database and return it
    names = list(pd.unique(vmh.loc[vmh["abbreviation"]==entry,"fullName"]))

    return(annotations, names)


def annotateVMH_id(metabolites: list[MeMoMetabolite]) -> int:
    """
    Annotate a list of metabolites with the entries from VMH. Look for VMH ids in the metabolite._id slot and if one is found use these.
    """

    # load the database
    config = get_config()
    db_path =  os.path.join(get_database_path(), config["databases"]["VMH"]["file"])
    with open(db_path) as f:
        vmh_json = json.load(f)
    vmh_json = json.dumps(vmh_json['results'])
    vmh = pd.read_json(StringIO(vmh_json))
    
    counter = 0
    for met in metabolites:
        new_met_anno = dict()
        #if met._id in bigg["universal_bigg_id"]: # does not work, for whatever reason
        if any(vmh["abbreviation"]==met._id):
            new_met_anno,new_names = annotateVMH_entry(entry = met._id,
                    annotations = new_met_anno,
                    database = vmh)

            # remove duplicates from the annotations
            for key in new_met_anno.keys():
               new_val = list(set(new_met_anno[key]))
               new_met_anno[key] = new_val
        
            # add the annotations to the slot in the metabolites
            if len(new_met_anno) > 0:
                met.add_annotations(new_met_anno)
                counter = counter +1
            
            # should count inchi?
            
            # add new names
            if len(new_names) > 0:
                met.add_names(new_names)

    # return a message how many annotations could be found
    print("Found " + str(counter) + " annotations in the VMH database, which have been added")

    return counter


def annotateVMH(metabolites: list[MeMoMetabolite]) -> int:
    """
    Annotate a list of metabolites with the entries from VMH. Look for VMH ids in the annotation slot of the metabolites and if one is found use these.
    """
    
    # load the database
    config = get_config()
    db_path =  os.path.join(get_database_path(), config["databases"]["VMH"]["file"])
    with open(db_path) as f:
        vmh_json = json.load(f)
    vmh_json = json.dumps(vmh_json['results'])
    vmh = pd.read_json(StringIO(vmh_json))
    
    counter = 0
    # go through the metabolites and check if there is data which can be added
    for met in metabolites:
        new_met_anno = dict()
        if "vmhmetabolite" in met.annotations.keys():
            for entry in met.annotations["vmhmetabolite"]:
                new_met_anno, new_names = annotateVMH_entry(entry = entry,
                        annotations = new_met_anno,
                        database = vmh)

            # remove duplicates from the annotations
            for key in new_met_anno.keys():
               new_val = list(set(new_met_anno[key]))
               new_met_anno[key] = new_val
            
            # add the annotations to the slot in the metabolites
            if len(new_met_anno) > 0:
                met.add_annotations(new_met_anno)
                counter = counter +1
            
            # should count inchi?
            
            # add new names
            if len(new_names) > 0:
                met.add_names(new_names)

    # return a message how many annotations could be found
    print("Found " + str(counter) + " annotations in the VMH database, which have been added")

    return counter

