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

    # inchis are not available but inchikeys are, that might take a while
    inchis = []
    for url in tqdm(dat.database_links.fillna(""), desc = "Translating inchi keys", unit = "item"):
        if url != "": 
            inchi_key = getInchiKeyFromURLS(url)
            if inchi_key == None:
                inchis.append("")
            else:
                inchi = inchikeyToInchi(inchi_key)
                if inchi != None:
                    inchis.append(inchi)
                else:
                    inchis.append("")
    inchis_series = pd.Series(inchis, index = dat.bigg_id, name = "inchi")



