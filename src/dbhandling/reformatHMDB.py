
# Enables execution as python  src/dbhandling/reformatHMDB.py Databases/hmdb_metabolites.zip hmdb_metabolites.xml
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

import sys

import pandas as pd
from src.dbhandling.reformatAux import *
import warnings

def concatCols(
    x: pd.Series,
    colNames: list[str],
    sep: str = '|'
) -> str:
    """
    Concatenate multiple fields from a pandas Series,
    ignoring empty strings and None values,
    joined by `sep`.
    
    Parameters:
        x (pd.Series): a row from a DataFrame
        fields (list[str]): list of column names to concatenate
        sep (str): separator string
    
    Returns:
        str: concatenated string of non-empty fields
    """
    parts = []
    for field in colNames:
      val = x.get(field)
      if val is not None and val != '':
          parts.append(str(val))

    for p in parts:
      original = p
      index: int=  p.find(sep)
      if index > 0:
        p = p.replace(sep,"*")
        print(f"Replacing {sep} in {original}")
        
    return sep.join(parts)

_keys: dict[str, str] = {
    "accession" : "id", 
    "name" : "name",
     "synonyms" : "synonyms",  # Empty column, might chang in the future
    "chemical_formula": "formula",
    "iupac_names" : "iupac_names",
    "traditional_iupac" : "traditional_iupac",
    "smiles" : "smiles",
    "inchi" : "inchi",
    "inchikey" : "inchi_key",
    "chemspider_id" : "chemspider",
    "drugbank_id" :  "drugbank",
    "metlin_id" :  "metlin",
    "foodb_id" : "food.compound",
    "pubchem_compound_id" : "pubchem.compound", 
    "chebi_id": "chebi", 
    "kegg_id" : "kegg.compound",
    "biocyc_id" : "biocyc", 
    "bigg_id": "bigg.metabolite", 
    "vmh_id" : "vmhmetabolite",
}
anno_keys: dict[str, str] = {
       "chemspider" : "chemspider",
       "drugbank" :  "drugbank",
       "metlin" :  "metlin",
       "food.compound" : "food.compound",
       "pubchem.compound" : "pubchem.compound",
       "chebi": "chebi",
       "kegg.compound" : "kegg.compound",
       "biocyc" : "biocyc",
       "bigg.metabolite": "bigg.metabolite",
       "vmhmetabolite" : "vmhmetabolite",
   }



def getAnnosPerEntry(dat:pd.DataFrame, met_id: str) -> dict[str,list[str]]:
    dat_sel = dat.loc[dat.id == met_id,]

    anno = {"hmdb":[met_id]}
    for key in anno_keys.keys():
      # Get whole key's col and get the cell (there is only one .... hopefully)
      val = dat_sel[key].iloc[0]
      if val != None and val != "":
          anno.setdefault(anno_keys[key], []).append(val)
    return(anno)

def getAnnos(dat:pd.DataFrame) -> pd.Series:
    """takes the data frame frwom VMH request and returns a pandas series of strings. Each string is a dictionary for the database annotations and can be simply evaluated (eval()) to be transformed in the correct format"""
    anno_series = pd.Series()
    for i in range(len(dat)):
        print(i, len(dat))
        met_id = dat.iloc[i]["id"]
        anno = getAnnosPerEntry(dat, met_id)
        anno_series[len(anno_series)] = str(anno)


    anno_series.index = list(dat["id"])
    return(anno_series)

def rename_columns_safe(df, rename_dict):
  for old_name, new_name in rename_dict.items():
    if old_name in df.columns:
      df.rename(columns={old_name: new_name}, inplace=True)
    else:
      print(f"Warning: Column '{old_name}' not found in DataFrame — skipping rename.")
  return df


def main():
    df = getData("HMDB")

    columns_to_concat = ["name", "synonyms", "iupac_names", "traditional_iupac"]
    # Apply the function row-wise to create a new column
    df['name'] = df.apply(concatCols, axis=1, colNames=columns_to_concat, sep='_|_')
    print("Added new col")
    if df is None:
        warnings.warn("Error during concatenation of columns in HMDB. Not Database created.")
        return 1
    df.drop(columns = ["synonyms", "iupac_names", "traditional_iupac"])
    print("droped old cols")
    # Don't ask me why, just df.rename(_keys) did not work
    rename_columns_safe(df, _keys)
    print("get annos")
    DBs = getAnnos(df)
    dat_all = df.loc[:,["id","name","inchi"]]
    dat_all.index = dat_all["id"]
    DBs = DBs.astype(str)
    DBsF = DBs.to_frame()
    print("starting join")
    dat_all = pd.concat([dat_all, DBsF], axis = 1)
    dat_all.columns = ["id", "name","inchi","DBs"]
    df.index = df["id"]
    dat_all = pd.concat([dat_all, df[["smiles", "inchi_key"]]], axis = 1)
    writeData(dat_all, db = "HMDB")


if __name__ == '__main__':
    # Run as python -m src.dbhandling.reformatHMDB 
    # or as python  src/dbhandling/reformatHMDB.py 
    main()

