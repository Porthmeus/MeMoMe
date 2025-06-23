# Porthmeus
# 22.05.25
from src.download_db import get_config, get_database_path
import pandas as pd
import os
import json

from io import TextIOWrapper
import zipfile
import xml.etree.ElementTree as ET
import csv


def strip_namespace(tag: str ) -> str:
    # Remove namespace from tag: {namespace}tag -> tag
    # Everything will have  http://www.hmdb.ca} preprended
    if '}' in tag:
        return tag.split('}', 1)[1]
    return tag

def xml_to_pandas_lazy(xml_stream: TextIOWrapper, record_tag: str):
    # TODO Merge name and synonyms and iupac_names
    # Translate to identifiers.org prefix
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


    i = 0
    rows: list[dict[str, list[str] | str]] = []
    for _, elem in ET.iterparse(xml_stream, events=('end',)):
        i += 1
        tag = strip_namespace(elem.tag)
        if tag == record_tag:
            row = {strip_namespace(child.tag): (child.text or '').strip() for child in elem}
            filtered_row: dict[str, list[str] | str] = {k: row.get(k, '') for k in _keys.keys()}
            rows.append(filtered_row)
            elem.clear()  # Free memory

    # Convert to list to staisfy type checker 
    df = pd.DataFrame(rows, columns=list(_keys.keys()))
    return df

def __json_to_dataframe(path) -> pd.DataFrame:
  with open(path, "r") as f: 
    vmh = json.load(f)["results"]
    return(pd.DataFrame(vmh))

def getData(db:str) -> pd.DataFrame:
    # loads database and the indentifier prefixes and returns them
    # get the database
    config = get_config()
    # no HMDB, can't be used because we do not want to load the whole DB
    dbs_csv = ["BiGG","ModelSeed"]
    dbs_json = ["VMH"]
    dbs_xml = ["HMDB"]
    path = os.path.join(get_database_path(),config["databases"][db]["file"])
    if db in dbs_csv:
        dat = pd.read_csv(path,
                          sep = "\t",
                          low_memory=False)
    elif db in dbs_json:
        file = os.path.join(path)
        dat = __json_to_dataframe(file)
    elif db in dbs_xml:
      with zipfile.ZipFile(path, 'r') as z:
        print(z.namelist())
        with z.open("hmdb_metabolites.xml") as xml_file:
            xml_text = TextIOWrapper(xml_file, encoding='utf-8')
            dat = xml_to_pandas_lazy(xml_text, "metabolite")
    else:
        raise ValueError("db must be one of " + str(dbs_csv + dbs_json))
    return(dat)


def writeData(dat:pd.DataFrame, db:str) -> None:
    config = get_config()
    outfile = os.path.join(get_database_path(), config["databases"][db]["reformat"])
    dat.to_csv(outfile)

