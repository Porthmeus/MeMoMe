import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

import sys

import pandas as pd
from src.dbhandling.reformatAux import *
import warnings



def main():
    df = getData("ChEBI")
    df.index = df["CHEBI_ID"]
    df = df.rename(columns = {'CHEBI_ID': "id", "InChI": "inchi"})

    df.index.rename('', inplace=True)
    return df


if __name__ == '__main__':
    # Run as python -m src.dbhandling.reformatHMDB 
    # or as python  src/dbhandling/reformatHMDB.py 
    df = main()
    writeData(df, "ChEBI")
