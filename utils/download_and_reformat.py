#TODO use configf.yml location
#TODO USE Path
import pandas as pd
import json
import os, sys

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)

db_location = "Databases/"
from src.download_db import *
from src.annotation.annotateAux import load_database

def __json_to_tsv(path: str) -> pd.DataFrame:
  with open(path, "r") as f: 
    vmh = json.load(f)["results"]
    return(pd.DataFrame(vmh))


def __load_csv(path: str) -> pd.DataFrame:
  df = pd.read_csv(path, low_memory=False)
  return(df)


def __load_tsv(path: str) -> pd.DataFrame:
  df = pd.read_csv(path, low_memory=False, sep = "\t")
  return(df)

#vmh = load_database(get_config()["databases"]["VMH"]["file"], allow_missing_dbs = False, conversion_method= __json_to_tsv)
## Relevant VMH columns
#keys = ["abbreviation", "fullName", "biggId", "keggId", "cheBlId", "inchiString", "inchiKey", "smile", "hmdb", "metanetx", "seed", "biocyc"]
#vmh = vmh[keys]
#vmh.to_csv(db_location + "vmh_formatted.csv")
#
#hmdb = load_database(get_config()["databases"]["HMDB"]["file"], allow_missing_dbs = False, conversion_method= __load_csv)
## Relevant HMDB columns
#keys = ["accession",  "chemical_formula", "iupac_name", "traditional_iupac", "smiles", "inchi", "chebi_id", "pubchem_compound_id", "kegg_id", "biocyc_id", "bigg_id", "vmh_id"]
#hmdb = hmdb[keys]
#hmdb.to_csv(db_location + "hmdb_formatted.csv")


#bigg = load_database(get_config()["databases"]["BiGG"]["file"], allow_missing_dbs = False, conversion_method= __load_tsv)
#keys = ["universal_bigg_id", "database_links"]
#bigg = bigg[keys]
#bigg.to_csv(db_location + "bigg_formatted.csv")
#
#
#chebi = load_database(get_config()["databases"]["ChEBI"]["file"], allow_missing_dbs = False, conversion_method= __load_tsv)
#keys = ["CHEBI_ID", "InChI"]
#chebi = chebi[keys]
#chebi.to_csv(db_location + "chebi_formatted.csv")
#print(chebi.columns)


seed = load_database(get_config()["databases"]["ModelSeed"]["file"], allow_missing_dbs = False, conversion_method= __load_tsv)
print(seed.columns)
