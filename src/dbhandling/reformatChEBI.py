import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

import sys

import pandas as pd
from src.dbhandling.reformatAux import *
import warnings



def reformatChEBI(df: pd.DataFrame):
    # Prepend "CHEBI:" to the CHEBI_ID column
    df["CHEBI_ID"] = "CHEBI:" + df["CHEBI_ID"].astype(str)

    # Set as index
    df.index = df["CHEBI_ID"]

    # Rename columns
    df = df.rename(columns={'CHEBI_ID': "id", "InChI": "inchi"})

    # Remove index name
    df.index.rename('', inplace=True)

    return df   


if __name__ == '__main__':
    # Run as python -m src.dbhandling.reformatHMDB 
    # or as python  src/dbhandling/reformatHMDB.py 
    df = getData("ChEBI")
    df_ref = reformatChEBI(df)
    writeData(df_ref, "ChEBI")
