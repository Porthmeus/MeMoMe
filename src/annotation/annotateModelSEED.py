# Porthmeus | None
# 07.07.23

'''
This file contains function to annotated MeMoMetabolites in a bulk manner. This means a list of MeMoMetabolites is parsed and the annotation takes place on the metabolites in that list. This is mainly to use within the MeMoModel.annotate() function
'''

import os
import pandas as pd
import json
import warnings
from io import StringIO
import re
from src.annotation.annotateInchiRoutines import findOptimalInchi, smile2inchi
from src.MeMoMetabolite import MeMoMetabolite
from src.download_db import get_config, get_database_path
from src.annotation.annotateAux import AnnotationResult, load_database, handleIDs, handleMetabolites
import ast
import re

# load the prefixes for the identifies.org list that is needed for the correction of the key in the annotation dictionary
config = get_config()
identifiers_file = os.path.join(get_database_path(), config["databases"]["Identifiers"]["file"])


def extractModelSEEDAnnotationsFromAlias(alias:str):
  '''
  Takes an entry from the model seed database ("https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/compounds.tsv") and retrieves the crosslinked databases into the annotation dictionary
  '''
  alias = alias.replace('"',"''")
  # Each type of identifier is separated by | i.e. we have Name: .... | BiGG: ... |
  # Now each list entry is corresponds to one db 
  ret1 = alias.split("|")

  result = {}
  # A entry looks like this DB: name1; name2; .... name_n 
  for entry in ret1:
    # Thus we get the part before the : and use it as a dict key
    key, value = entry.split(':')  
    # The part after the : contains the spafe to the right from the : and the ;.
    # Thus we remov all " " aka whitespace 
    # Thus name1; name2; .... name_n turns into name1;name2;name3; ...;name_n
    # Then split by ; to get all the names
    # This line removes the whitespace after the :
    value = value.strip()
    # This line removes the whitespace after the ;
    result[key.strip()] = re.sub(r";\s", ";", value).split(";")

  names = result.pop("Name", [])
  return result, names

#def correctAnnotationKeys(anno:dict, allow_missing_dbs: bool = False) -> dict:
#    "Takes an annotation dictionary and tries to correct the keys in the dictionary to conform with the identifiers.org prefixes"
#
#    ## TODO: currently the data base is read everytime we run this function (which could be several k times) - find a better solution for that
#    # try to find the correct identifiers from identifiers.org
#    ## first get possible prefixes
#    prefixes = identifier_prefixes
#
#    ## now check if the keys of the annos can be found in the prefixes
#    new_anno_corrected = dict()
#    for key, value in anno.items():
#        if key.lower() in prefixes:
#            new_anno_corrected[key.lower()] = value
#        # handle kegg drug vs. compound
#        elif key.lower() == "kegg":
#            # sort the values into compounds and drugs
#            drugs = []
#            compounds = []
#            for val in value:
#                if val.startswith("D"):
#                    drugs.append[val]
#                elif val.startswith("C"):
#                    compounds.append[val]
#            if len(drugs) >0:
#                new_anno_corrected["kegg.drug"] = drugs
#            if len(compounds) >0:
#                new_anno_corrected["kegg.compound"] = compounds
#        # remove undefined pubchem entries - compounds and substances can have the same ID, but are usually different metabolites 
#        elif key.lower() == "pubchem":
#            warnings.warn("Removing pubchem entries because of missing information about compound or substance:" + str(value))
#        # handle database specific stuff like the .compound suffix
#        elif key.lower()+".compound" in prefixes:
#            new_anno_corrected[key.lower() + ".compound"] = value
#        elif key.lower()+".metabolite" in prefixes:
#            new_anno_corrected[key.lower() + ".metabolite"] = value
#        elif key.lower()+".cpd" in prefixes:
#            new_anno_corrected[key.lower() + ".cpd"] = value
#        elif key.lower()+".chemical" in prefixes:
#            new_anno_corrected[key.lower() + ".substance"] = value
#        elif key.lower()+".substance" in prefixes:
#            new_anno_corrected[key.lower() + ".chemical"] = value
#        elif key.lower()+".entity" in prefixes:
#            new_anno_corrected[key.lower() + ".entity"] = value
#        elif key.lower()+".ligand" in prefixes:
#            new_anno_corrected[key.lower() + ".ligand"] = value
#        elif key.lower()+".smallmolecule" in prefixes:
#            new_anno_corrected[key.lower() + ".ligand"] = value
#        elif key.lower()+".inhibitor" in prefixes:
#            new_anno_corrected[key.lower() + ".inhibitor"] = value
#        elif key.lower()+".pthcmp" in prefixes:
#            new_anno_corrected[key.lower() + ".pthcmp"] = value
#        
#    return(new_anno_corrected)

def handle_seed_entry(aliases: pd.DataFrame) -> tuple[dict, list, str]:
# get the information from the database
  new_anno = dict()
  new_names = list()
  if len(aliases) >0:
      for alias in aliases:
          annos, names = extractModelSEEDAnnotationsFromAlias(alias)
          # combine annos for the entry
          for key,value in annos.items():
              if key in new_anno.keys():
                  new_anno[key].extend(value)
              else:
                  new_anno[key] = value
          # combine names
          new_names.extend(names)
      
      # remove duplicates
      for key in new_anno.keys():
          new_anno[key] = list(set(new_anno[key]))
      new_names = list(set(new_names))

  return(new_anno, new_names, "Seed")


def annotateModelSEED_entry(entry:str,  database:pd.DataFrame = pd.DataFrame(), allow_missing_dbs: bool = False) -> tuple[dict, list, str]:
    """
    A small helper function to avoid redundant code
    Uses a ModelSEED identifiers and annotates it with the identifiers.org entries.
    The results are stored in the annotations dictionary. This is empty by
    default, but if it exists, the entries will be simply added (duplicates are
    allowed here). Further it will return a list of alternative names, and pKa/pKb values
    database - is either empty (default) or a pandas data.frame with the data from the modelseed homepage for the metabolites. If it is empty, the function will try to load the database from the config file
    """
    
    if len(database) == 0:
      mseed =  load_database(get_config()["databases"]["ModelSeed"]["file"], 
                            allow_missing_dbs, 
                            lambda path: pd.read_csv(path, sep="\t", low_memory=False))
    else:
      mseed = database

    if mseed.empty:
      return dict(), list(), "Seed"

    # extract the relevant annotation information and return it
    aliases = mseed.loc[mseed["id"]==entry,"aliases"]
    # check if there are entries which are not NA
    aliases = aliases.loc[~pd.isna(aliases)]
    return handle_seed_entry(aliases)
    


def annotateModelSEED_id(metabolites: list[MeMoMetabolite], allow_missing_dbs: bool = False) ->AnnotationResult:
    """
    Annotate a list of metabolites with the entries from ModelSeed. Look for ModelSeed ids in the metabolite._id slot and if one is found use these. Since ModelSeed does provide any InChI strings, pKs and pkB, all these elements will be added to the MeMoMetabolites, if possible"""

    
    mseed =  load_database(get_config()["databases"]["ModelSeed"]["file"], 
                          allow_missing_dbs, 
                          lambda path: pd.read_csv(path, sep="\t", low_memory=False))
    if mseed.empty:
      return AnnotationResult(0, 0,0 )
    
    return handleIDs(mseed, metabolites, "id", annotateModelSEED_entry)

def annotateModelSEED(metabolites: list[MeMoMetabolite], allow_missing_dbs: bool = False) ->AnnotationResult:
  """
  Annotate a list of metabolites with the entries from ModelSeed. Look for ModelSeed ids in the metabolite.anotation slot and if one is found use these. Since ModelSeed does provide any InChI strings, pKs and pkB, all these elements will be added to the MeMoMetabolites, if possible
  """


  # load the database
  mseed =  load_database(get_config()["databases"]["ModelSeed"]["file"], 
                        allow_missing_dbs, 
                        lambda path: pd.read_csv(path, sep="\t", low_memory= False))
  if mseed.empty:
    return AnnotationResult(0, 0,0 )


  return handleMetabolites(mseed, metabolites, "seed.compound", annotateModelSEED_entry)
