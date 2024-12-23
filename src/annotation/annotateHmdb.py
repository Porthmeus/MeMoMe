import json
import pandas as pd
from src.download_db import get_config
from src.annotation.annotateAux import  load_database, AnnotationResult, handleMetabolites
from src.MeMoMetabolite import MeMoMetabolite
import warnings

def normalize_hmdb_id(hmdb_id: str):
  """
  params: A hmdb id
  VMH models contain keys like HMDB00972, but in HMDB the id is HMDB0000972. I.e. the digit part is 7 long not 5.
  """
  prefix, numeric_part = hmdb_id[:4], hmdb_id[4:]
  if prefix.upper() != "HMDB":
    warnings.warn(f"The id ({hmdb_id}) is not a HMDB id. We will not look for matches. We are taking the first for characters ({prefix}) and capitalize it ({prefix.upper()}). If those are not HMDB we consider this an error.")
    return ""
  normalized_id = f"{prefix}{numeric_part.zfill(7)}"
  return normalized_id

def handle_HMDB_entries(hmdb, entry: str):
  # VMH models contain keys like HMDB00972, but in HMDB the id is HMDB0000972. I.e. the digit part is 7 long not 5.
  entry = normalize_hmdb_id(entry)
  # The db col is called inchi, but it's inchi strings
  keys = ["accession",  "chemical_formula", "iupac_name", "traditional_iupac", "smiles", "inchi", "chebi_id", "pubchem_compound_id", "kegg_id", "biocyc_id", "bigg_id", "vmh_id"]
  hmdb = hmdb.loc[hmdb["accession"] == entry,]

  pd.set_option('future.no_silent_downcasting', True)
  hmdb = hmdb.replace(to_replace=r'\n', value='', regex=True)
  fullName = list(set(hmdb["name"]))
  hmdb = hmdb[keys]
  annotations = dict()
  if len(hmdb) > 1:
    raise Exception("More than one value for key: " + entry)

  for index, row in hmdb.iterrows():
    for key, value in row.items():
      if (not pd.isna(value)) and (len(str(value)) > 0):
        # If key does not exist, create it and append the value
        annotations.setdefault(key, []).append(value)
  return(annotations, fullName)

def __load_csv(path) -> pd.DataFrame:
  df = pd.read_csv(path, low_memory=False)
  return(df)


def annotateHMDB_entry(entry: str,  database: pd.DataFrame = pd.DataFrame(), allow_missing_dbs: bool = False) -> tuple[dict, list]:
    # check if the database was given, if not, try to load it
    if len(database) > 0:
      hmdb = database
    else:
      hmdb = load_database(get_config()["databases"]["HMDB"]["file"], allow_missing_dbs, __load_csv)
    
    if hmdb.empty:
      return dict(), list()

    annotations, names = handle_HMDB_entries(hmdb, entry)
    return annotations, names



def annotateHMDB(metabolites: list[MeMoMetabolite], allow_missing_dbs: bool = False) -> AnnotationResult:
    """
    Annotate a list of metabolites with the entries from VMH. Look for VMH ids in the annotation slot of the metabolites and if one is found use these.
    The function will directly add the annotation and names to the MeMoMetabolite object.
    Return value: a tuple of 3 denoting the number of changes metabolites for the following values: [inchi_strings, annotations, names]
    """
    
    hmdb = load_database(get_config()["databases"]["HMDB"]["file"], allow_missing_dbs, __load_csv)
    
    if hmdb.empty:
      return AnnotationResult(0, 0, 0)

    return handleMetabolites(hmdb, metabolites, "hmdb", annotateHMDB_entry)

