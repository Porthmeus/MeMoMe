import pandas as pd
from src.download_db import get_config, get_database_path
from src.MeMoMetabolite import MeMoMetabolite
from src.annotation.annotateAux import AnnotationResult, load_database, handleIDs, handleMetabolites, EntryAnnotationFunction, DBName, AnnotationKey, DBKey
from typing import Callable, List
import json


# TODO What exactly is STUFF
def annotateDB(metabolites: list[MeMoMetabolite], db_name: DBName, annotation_key_name: DBKey, loading_function,  
               entry_annotation_function: EntryAnnotationFunction , allow_missing_dbs: bool = False) -> AnnotationResult:
  """
  Annotate a list of metabolites with the entries from DB. Look for DB ids in the annotation slot of the metabolites and if one is found use these. 
  The function will directly add the annotation and names to the MeMoMetabolite object.
  Return value: a tuple of 3 denoting the number of changes metabolites for the following values: [inchi_strings, annotations, names]
  """
  db = load_database(get_config()["databases"][db_name]["file"], allow_missing_dbs, loading_function)

  if db.empty:
    return AnnotationResult(0, 0, 0)
  return handleMetabolites(db, metabolites, annotation_key_name, entry_annotation_function)


def annotateDB_entry(entry: AnnotationKey, db_name: DBName, loading_function,  
                     handle_function: Callable[[pd.DataFrame, AnnotationKey], tuple[dict, list, str]],
                     database: pd.DataFrame = pd.DataFrame(), allow_missing_dbs: bool = False) -> tuple[dict, list, str]:
    """
    A small helper function to avoid redundant code
    Uses a DB identifiers and annotates it with the identifiers.org entries.
    database - is either empty (default) or a pandas data.frame with the data from the DB for the metabolites. If it is empty, the function will try to load the database from the config file
    Returns a tuple containing a dictionary for the extracted annotations and a list of new names for the metabolite
    """
    # check if the database was given, if not, try to load it
    if len(database) > 0:
      db = database
    else:
      db = load_database(get_config()["databases"][db_name]["file"], allow_missing_dbs, loading_function)
    
    if db.empty:
      return dict(), list(), db_name 

    annotations, names, source = handle_function(db, entry)
    return annotations, names, source


def annotateDB_id(metabolites: list[MeMoMetabolite], db_name: DBName, db_key: DBKey, loading_function: Callable[[str], pd.DataFrame], 
                   entry_annotation_function: EntryAnnotationFunction, allow_missing_dbs: bool = False) -> AnnotationResult:
  """
  Annotate a list of metabolites with the entries from DB . Look for DB ids in the metabolite._id slot and if one is found use these. 
  """
  # check if the database was given, if not, try to load it
  db = load_database(get_config()["databases"][db_name]["file"], allow_missing_dbs, loading_function)
  
  if db.empty:
    return AnnotationResult(0, 0, 0)
  

  return handleIDs(db, metabolites, db_key, entry_annotation_function)
