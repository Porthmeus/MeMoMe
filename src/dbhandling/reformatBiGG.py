# Porthmeus
# 20.01.25

# functions to reformat the BiGG database to a standard format which can be easily read  
from src.download_db import get_config, get_database_path
from src.parseMetaboliteInfos import getAnnoFromIdentifierURL
from src.annotation.annotateInchiRoutines import *
from tqdm import tqdm
from io import StringIO
from math import isnan
import re
import pandas as pd
import os
import json
import tempfile
import requests
import subprocess

def getData() -> pd.DataFrame:
    # loads the modelSeed database and the indentifier prefixes and returns them
    # get the database
    config = get_config()
    dat = pd.read_csv(os.path.join(get_database_path(),config["databases"]["BiGG"]["file"]),
                      sep = "\t",
                      low_memory=False)
    return(dat)


def handleDBs(urls:str) -> str:
    # reformat the list of URLs given by BiGG to a dictionary of database name (key) and database id (values)
    if urls !="":
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
    return(str(anno_dbs))

def getInchiKeyFromURLS(urls:str) -> str:
    # get the inchikey from a string if it exists
    pat = r"identifiers\.org/inchikey/([a-zA-Z0-9\-]+)"
    match = re.search(pat, urls)
    if match:
        inchi_key = match.group(1)
    else:
        inchi_key = None
    return(inchi_key)

def downloadPBCInchiKey2String() -> str:
    # create a temporary directory, download the database and return the path to the file
    url = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-InChI-Key.gz"
    with tempfile.NamedTemporaryFile(delete=False, mode='wb') as temp_file:
        with requests.get(url, stream=True) as response:
            response.raise_for_status()  # Ensure request was successful
            for chunk in tqdm(response.iter_content(chunk_size=8192), desc = "Downloading Pubchem Database"):  # Read in chunks of 8 KB
                temp_file.write(chunk)  # Write chunk to file
        
        temp_file_path = temp_file.name  # Get the file path
    return(temp_file_path)


def createInchiKey2String(pbc_table:str, inchiKeys:list[str]) -> str:
    # use zgrep to create a conversion table from inchi key to inchi string via the pubchem database


    # safe the inchi keys in a temporary file
    with tempfile.NamedTemporaryFile(delete=False, mode = "w") as temp_file:
        for entry in inchiKeys:
            if entry !="":
                temp_file.write(entry + "\n")
        keys_file = temp_file.name
    
    # create subprocess to zgrep the inchi keys
    result = subprocess.run(
        ["zgrep", "-f", keys_file, pbc_table], 
        text=True,  # Ensures output is returned as a string
        capture_output=True,  # Captures stdout and stderr
        check=True  # Raises an error if the command fails
    )
    
    return(result.stdout)  # Print matched lines
    os.remove(keys_file)

def writeData(dat:pd.DataFrame) -> None:
    config = get_config()
    outfile = os.path.join(get_database_path(), config["databases"]["BiGG"]["reformat"])
    dat.to_csv(outfile)

def reformatBiGG() -> None:
    # get names, DBs and Inchi columns
    dat = getData()

    # names
    names = dat.name.fillna(value = "").apply(lambda x: str(x).replace(";",","))
    names.index = dat.bigg_id
    names.name = "names"

    # dbs
    dbs = dat.database_links.fillna(value = "").apply(handleDBs)
    dbs.index = dat.bigg_id
    dbs.name = "DBs"

    # inchis are not available but inchikeys are, get these and translate them
    inchis_keys = []
    for url in dat.database_links.fillna(""):
        if url != "": 
            inchi_key = getInchiKeyFromURLS(url)
            if inchi_key == None:
                inchis_keys.append("")
            else:
                inchis_keys.append(inchi_key)
        else:
            inchis_keys.append("")

    # create a unique 
    inchi_keys = list(set(inchis_keys))
    # create a series
    inchis_keys =pd.DataFrame({"inchi_key":inchis_keys,
                               "bigg_id" : dat.bigg_id})
    
    # download the pubchem database and get the conversion table
    pbc_table = downloadPBCInchiKey2String()
    convert_table = createInchiKey2String(pbc_table, inchi_keys)
    tbl_conversion = pd.read_csv(StringIO(convert_table), delimiter = "\t",usecols=[1,2], header =None)
    tbl_conversion.columns = ["inchi","inchi_key"]
    # create a proper table and r
    inchis = pd.merge(inchis_keys, tbl_conversion,on ="inchi_key", how = "left")
    inchis.index = inchis.bigg_id
    # 
    dat_all = pd.concat([names, dbs], axis = 1)
    dat_all = pd.merge(dat_all, inchis, left_index =True, right_index = True)
    
    # write the database
    writeData(inchi_keys)

    # remove tempfiles
    os.remove(pbc_table)
