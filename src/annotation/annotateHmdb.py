import json
import pandas as pd
from src.download_db import get_config
from src.annotation.annotateAux import  load_database, AnnotationResult, handleMetabolites, DBName, DBKey
from src.MeMoMetabolite import MeMoMetabolite
import warnings
from src.annotation.annotateALL import *



def normalize_hmdb_id(hmdb_id: AnnotationKey) -> AnnotationKey:
  """
  params: A hmdb id
  VMH models contain keys like HMDB00972, but in HMDB the id is HMDB0000972. I.e. the digit part is 7 long not 5.
  """
  prefix, numeric_part = hmdb_id[:4], hmdb_id[4:]
  if prefix.upper() != "HMDB":
    warnings.warn(f"The id ({hmdb_id}) is not a HMDB id. We will not look for matches. We are taking the first for characters ({prefix}) and capitalize it ({prefix.upper()}). If those are not HMDB we consider this an error.")
    return AnnotationKey("")
  normalized_id = f"{prefix}{numeric_part.zfill(7)}"
  return AnnotationKey(normalized_id)

def handle_HMDB_entries(hmdb: pd.DataFrame, entry: AnnotationKey):

  # VMH models contain keys like HMDB00972, but in HMDB the id is HMDB0000972. I.e. the digit part is 7 long not 5.
  entry = normalize_hmdb_id(entry)
  # The db col is called inchi, but it's inchi strings
  keys = ["accession",  "chemical_formula", "iupac_name", "traditional_iupac", "smiles", "inchi", "chebi_id", "pubchem_compound_id", "kegg_id", "biocyc_id", "bigg_id", "vmh_id"]
  hmdb = hmdb.loc[hmdb["accession"] == entry,]
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

def __load_csv(path: str) -> pd.DataFrame:
  df = pd.read_csv(path, low_memory=False)
  return(df)


def annotateHMDB_entry(entry: AnnotationKey,  database: pd.DataFrame = pd.DataFrame(), allow_missing_dbs: bool = False) -> tuple[dict, list, str]:
  return annotateDB_entry(entry, DBName("HMDB"), __load_csv, handle_HMDB_entries, database, allow_missing_dbs)

def annotateHMDB(metabolites: list[MeMoMetabolite], allow_missing_dbs: bool = False) -> AnnotationResult:
  return  annotateDB(metabolites, DBName("HMDB"), DBKey("HMDB"), __load_csv, annotateHMDB_entry, allow_missing_dbs)
