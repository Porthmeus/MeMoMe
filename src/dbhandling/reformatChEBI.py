import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
import sys
import pandas as pd
from src.dbhandling.reformatAux import *
import warnings
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def inchi_to_formula(inchi):
    try:
        mol = Chem.MolFromInchi(inchi)
        if mol is None:
            return None
        return rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        return None


def reformatChEBI(df: pd.DataFrame):
    # Create formula column from InChI
    df["formula"] = df["InChI"].apply(inchi_to_formula)

    # Prepend "CHEBI:" to the CHEBI_ID column
    df["CHEBI_ID"] = "CHEBI:" + df["CHEBI_ID"].astype(str)

    # Set as index
    df.index = df["CHEBI_ID"]

    # Rename columns
    df = df.rename(columns={
        "CHEBI_ID": "id",
        "InChI": "inchi"
    })

    # Remove index name
    df.index.rename('', inplace=True)

    return df



if __name__ == '__main__':
    # Run as python -m src.dbhandling.reformatHMDB 
    # or as python  src/dbhandling/reformatHMDB.py 
    df = getData("ChEBI")
    df_ref = reformatChEBI(df)
    writeData(df_ref, "ChEBI")
