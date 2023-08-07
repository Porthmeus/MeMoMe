import gzip
import sys

from rdkit import Chem
from rdkit.Chem import inchi
import time

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
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
    start = time.time()

    gzipped_sdf_file_path = sys.argv[1]
    print(f"Loading {gzipped_sdf_file_path}")
    molecules = read_gzipped_sdf_file(gzipped_sdf_file_path)

    # Now you can work with the molecules, for example, print the number of molecules:
    print("Number of molecules in gzipped SDF file:", len(molecules))
    id_to_inchi = {}
    for i in molecules:
        name = i.GetProp("PUBCHEM_SUBSTANCE_ID")  # Replace "PUBCHEM_COMPOUND_NAME" with the desired property name.
        inchi_str = inchi.MolToInchi(i)
        id_to_inchi[name] = inchi_str
    end = time.time()
    print(end - start)
    # You can also access properties of each molecule as shown in the previous example.
    print("TESt")
