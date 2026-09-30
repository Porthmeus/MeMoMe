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


import sys
import os
import warnings
import pandas as pd
from tqdm import tqdm

if __name__ == "__main__":
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from src.dbhandling.reformatAux import getData, writeData
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import RDLogger

# 1. Silence RDKit C++ warnings/errors completely
RDLogger.DisableLog('rdApp.*')

def inchi_to_formula(inchi):
    if pd.isna(inchi) or not isinstance(inchi, str) or inchi == "":
        return None
    try:
        mol = Chem.MolFromInchi(inchi)
        if mol is None:
            return None
        return rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        return None

def reformatChEBI(df: pd.DataFrame):
    # 2. Register tqdm with pandas to see a live progress bar
    tqdm.pandas(desc="Calculating ChEBI formulas")
    
    # Create formula column from InChI safely
    df["formula"] = df["InChI"].progress_apply(inchi_to_formula)

    # Prepend "CHEBI:" to the CHEBI_ID column
    df["CHEBI_ID"] = df["CHEBI_ID"].astype(str)

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
    df = getData("ChEBI")
    df_ref = reformatChEBI(df)
    writeData(df_ref, "ChEBI")
