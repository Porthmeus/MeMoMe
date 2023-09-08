print("test")
print("Current Python interpreter:")


import gzip
import sys
print(sys.executable)
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import inchi
import time
import pandas as pd
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
print("Imports are done")
def read_gzipped_sdf_file(file_path):
    molecules = []
    with gzip.open(file_path, 'rt') as f:
        # Read the decompressed content into a string
        sdf_content = f.read()
        # Use SDMolSupplier to read the SDF data from the string
        suppl = Chem.SDMolSupplier()
        suppl.SetData(sdf_content)
        for mol in suppl:
            if mol is not None:
                molecules.append(mol)
    return molecules


if __name__ == "__main__":
    print("args" + str(sys.argv))
    start = time.time()
    if len(sys.argv) <= 2:
        print("Not enough args supplied")
        sys.exit(1)
    gzipped_sdf_file_path = sys.argv[1]
    print(f"Loading {gzipped_sdf_file_path}")
    molecules = read_gzipped_sdf_file(gzipped_sdf_file_path)

    # Now you can work with the molecules, for example, print the number of molecules:
    print("Number of molecules in gzipped SDF file:", len(molecules))
    id_to_inchi = {}
    for i in molecules:
        name = i.GetProp("PUBCHEM_SUBSTANCE_ID")  # Replace "PUBCHEM_COMPOUND_NAME" with the desired property name.
        try:
            inchi_str = inchi.MolToInchi(i)
        except:
            inchi_str = "NA"
        finally:
            id_to_inchi[name] = inchi_str

    end = time.time()
    df = pd.DataFrame.from_dict(id_to_inchi, orient='index', columns=['inchi'])
    out = sys.argv[2] 
    print(out)
    df.to_csv(out)
    print(f"Took {end - start} seconds")
    # You can also access properties of each molecule as shown in the previous example.
