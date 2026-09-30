# Porthmeus
# 20.01.25

# functions to reformat the BiGG database to a standard format which can be easily read

import sys
import os

if __name__ == "__main__":
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from src.parseMetaboliteInfos import getAnnoFromIdentifierURL
from src.annotation.annotateInchiRoutines import *
from src.dbhandling.reformatAux import getData, writeData

from tqdm import tqdm
from io import StringIO
from math import isnan
import re
import pandas as pd
import tempfile
import requests
import subprocess
import warnings

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import RDLogger

# Silence noisy RDKit C++ backend warnings/errors
RDLogger.DisableLog('rdApp.*')


def inchi_to_formula(inchi):
    """
    Convert InChI string to molecular formula using RDKit.
    """
    if pd.isna(inchi) or not isinstance(inchi, str) or inchi == "":
        return None

    try:
        mol = Chem.MolFromInchi(inchi)

        if mol is None:
            return None

        return rdMolDescriptors.CalcMolFormula(mol)

    except Exception:
        return None


def handleDBs(urls: str) -> str:
    # reformat the list of URLs given by BiGG to a dictionary of database name (key) and database id (values)

    if urls != "":
        urls = urls.split(";")
        anno_dbs = dict()

        for url in urls:
            split_url = getAnnoFromIdentifierURL(url)

            if split_url[0] in anno_dbs.keys():

                # do not add duplicates
                if split_url[1] not in anno_dbs[split_url[0]]:
                    anno_dbs[split_url[0]].append(split_url[1])
                    anno_dbs[split_url[0]].sort()

            else:
                anno_dbs[split_url[0]] = [split_url[1]]

    else:
        anno_dbs = dict()

    return str(anno_dbs)


def getInchiKeyFromURLS(urls: str) -> str:
    # get the inchikey from a string if it exists

    pat = r"identifiers\.org/inchikey/([a-zA-Z0-9\-]+)"
    match = re.search(pat, urls)

    if match:
        inchi_key = match.group(1)
    else:
        inchi_key = None

    return inchi_key


def downloadPBCInchiKey2String(
    url: str = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-InChI-Key.gz"
) -> str:

    # create a temporary directory, download the database and return the path to the file

    with tempfile.NamedTemporaryFile(delete=False, mode='wb') as temp_file:

        with requests.get(url, stream=True) as response:

            response.raise_for_status()

            for chunk in tqdm(
                response.iter_content(chunk_size=1024 * 1024),
                desc="Downloading Pubchem Database (~6800it)"
            ):
                temp_file.write(chunk)

        temp_file_path = temp_file.name

    return temp_file_path


def createInchiKey2String(pbc_table: str, inchiKeys: list[str]) -> str:
    """
    Use zgrep to create a conversion table from InChIKey to InChI.
    """

    with tempfile.NamedTemporaryFile(delete=False, mode="w") as temp_file:

        for entry in inchiKeys:
            if entry != "":
                temp_file.write(entry + "\n")

        keys_file = temp_file.name

    result = subprocess.run(
        ["zgrep", "-f", keys_file, pbc_table],
        text=True,
        capture_output=True,
        check=True
    )

    os.remove(keys_file)

    return result.stdout


def joinValues(x: list[str], sep: str = "_|_") -> str:
    """
    Helper function to remove nan from string columns.
    """

    x = [val for val in x if val is not None and not pd.isna(val) and val != ""]
    ret = "_|_".join(sorted(set(x)))

    return ret


def sortDataFrame(bigg_data: pd.DataFrame) -> pd.DataFrame:
    """
    Removes duplicates caused by metabolites occurring in multiple compartments.
    """

    bigg_data = bigg_data.loc[
        :,
        ["universal_bigg_id", "name", "database_links", "old_bigg_ids"]
    ]

    bigg_data = bigg_data.loc[~bigg_data.duplicated(), :]

    bigg_combined = bigg_data.groupby("universal_bigg_id").agg(
        lambda x: joinValues(x, sep="_|_")
    )

    bigg_combined["old_bigg_ids"] = (
        bigg_combined["old_bigg_ids"]
        .apply(lambda x: "; ".join(sorted(set(x.replace("_|_", "; ").split("; ")))))
    )

    bigg_combined["database_links"] = (
        bigg_combined["database_links"]
        .apply(lambda x: "; ".join(sorted(set(x.replace("_|_", "; ").split("; ")))))
    )

    bigg_combined["universal_bigg_id"] = bigg_combined.index

    return bigg_combined


def reformatBiGG(
    dat: pd.DataFrame,
    inchiKey2String_url: str = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-InChI-Key.gz"
) -> pd.DataFrame:
    """
    Reformats BiGG database and returns a dataframe with:
    id, names, inchis, formulas and DB annotations.
    """

    dat = sortDataFrame(dat)

    # Names
    names = (
        dat.name
        .fillna("")
        .apply(lambda x: str(x).replace(";", ","))
    )

    names.index = dat.universal_bigg_id
    names.name = "name"

    # DBs
    dbs = (
        dat.database_links
        .fillna("")
        .apply(handleDBs)
    )

    dbs.index = dat.universal_bigg_id
    dbs.name = "DBs"

    # Extract InChIKeys
    inchis_keys = []

    for url in dat.database_links.fillna(""):

        if url != "":
            inchi_key = getInchiKeyFromURLS(url)

            if inchi_key is None:
                inchis_keys.append("")
            else:
                inchis_keys.append(inchi_key)

        else:
            inchis_keys.append("")

    unique_inchi_keys = list(set(inchis_keys))

    inchis_keys = pd.DataFrame({
        "inchi_key": inchis_keys,
        "id": dat.universal_bigg_id
    })

    # Download PubChem mapping table
    pbc_table = downloadPBCInchiKey2String(
        url=inchiKey2String_url
    )

    convert_table = createInchiKey2String(
        pbc_table,
        unique_inchi_keys
    )

    tbl_conversion = pd.read_csv(
        StringIO(convert_table),
        delimiter="\t",
        usecols=[1, 2],
        header=None
    )

    tbl_conversion.columns = [
        "inchi",
        "inchi_key"
    ]

    # Merge InChIs
    inchis = pd.merge(
        inchis_keys,
        tbl_conversion,
        on="inchi_key",
        how="left"
    )

    inchis.index = inchis.id

    dat_all = pd.concat(
        [names, dbs],
        axis=1
    )

    dat_all = pd.merge(
        dat_all,
        inchis,
        left_index=True,
        right_index=True
    )

    # --------------------------------------------------
    # Calculate molecular formulas with a progress bar
    # --------------------------------------------------
    tqdm.pandas(desc="Calculating BiGG formulas")
    dat_all["formula"] = dat_all["inchi"].progress_apply(
        inchi_to_formula
    )
    # --------------------------------------------------

    os.remove(pbc_table)

    dat_all_start = dat_all.loc[
        :,
        ["id", "name", "inchi", "formula", "DBs"]
    ]

    dat_all_end = dat_all.drop(
        ["id", "name", "inchi", "formula", "DBs"],
        axis=1,
        errors="ignore"
    )

    dat_all = pd.concat(
        [dat_all_start, dat_all_end],
        axis=1
    )

    return dat_all


def main():

    if len(sys.argv) != 1:
        print("Usage: python reformatBigg.py")
        sys.exit(1)

    try:

        dat = getData("BiGG")

        refDB = reformatBiGG(dat)

        writeData(
            refDB,
            db="BiGG"
        )

        print(
            "BiGG reformatting completed successfully. "
            "Output written to standard location."
        )

    except Exception as e:

        warnings.warn(
            f"Error during BiGG reformatting: {e}"
        )

        sys.exit(1)


if __name__ == "__main__":
    main()
