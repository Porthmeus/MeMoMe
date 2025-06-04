# Enables execution as python  src/dbhandling/reformatHMDB.py Databases/hmdb_metabolites.zip hmdb_metabolites.xml
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

import sys
import zipfile
import xml.etree.ElementTree as ET
import csv
from io import TextIOWrapper
import pandas as pd
from src.dbhandling.reformatAux import *
import warnings

 # TODO Merge name and synonyms and iupac_names
 # Translate to identifiers.org prefix
_keys: dict[str, str] = {
       "accession" : "hmdb", 
       "name" : "name",
        "synonyms" : "synonyms",  # Empty column, might chang in the future
       "chemical_formula": "formula",
       "iupac_names" : "iupac_names", 
       "traditional_iupac" : "traditional_iupac", 
       "smiles" : " smiles",
       "inchi" : "inci",
       "inchi_key" : "inchi_key", 
       "chemspider" : "chemspider",
       "foodb_id" : "food.compound",
       "pubchem_compound_id" : "pubchem.compound", 
       "chebi_id": "chebi", 
       "kegg_id" : "kegg.compound",
       "biocyc_id" : "biocyc", 
       "bigg_id": "bigg.metabolite", 
       "vmh_id" : "vmhmetabolite",
   }

def concatCols(
    x: pd.Series,
    colNames: list[str],
    sep: str = '|'
) -> str:
    """
    Concatenate multiple fields from a pandas Series,
    ignoring empty strings and None values,
    joined by `sep`.
    
    Parameters:
        x (pd.Series): a row from a DataFrame
        fields (list[str]): list of column names to concatenate
        sep (str): separator string
    
    Returns:
        str: concatenated string of non-empty fields
    """
    parts = []
    for field in colNames:
      val = x.get(field)
      if val is not None and val != '':
          parts.append(str(val))

    for p in parts:
      index: int=  p.find(sep)
      if index > 0:
        print(index)
        p = p.replace(sep,"*")
        print(f"Replacing {sep} in {p}")
        
    return sep.join(parts)

def strip_namespace(tag: str ) -> str:
    # Remove namespace from tag: {namespace}tag -> tag
    # Everything will have  http://www.hmdb.ca} preprended
    if '}' in tag:
        return tag.split('}', 1)[1]
    return tag

def xml_to_pandas_lazy(xml_stream: TextIOWrapper, record_tag: str):

    rows: list[dict[str, list[str] | str]] = []
    for _, elem in ET.iterparse(xml_stream, events=('end',)):
        tag = strip_namespace(elem.tag)
        if tag == record_tag:
            row = {strip_namespace(child.tag): (child.text or '').strip() for child in elem}
            filtered_row: dict[str, list[str] | str] = {k: row.get(k, '') for k in _keys.keys()}
            rows.append(filtered_row)
            elem.clear()  # Free memory

    # Convert to list to staisfy type checker 
    df = pd.DataFrame(rows, columns=list(_keys.keys()))
    return df

def main():
    if len(sys.argv) != 3:
        print("Usage: python lazy_xml_to_csv.py <zipfile.zip> <xml_filename_inside_zip>)")
        sys.exit(1)
    # Cant use get data because we do not want to read the whole file
    zip_path: str = sys.argv[1]
    xml_name: str = sys.argv[2]
    with zipfile.ZipFile(zip_path, 'r') as z:
        with z.open(xml_name) as xml_file:
            xml_text = TextIOWrapper(xml_file, encoding='utf-8')
            df = xml_to_pandas_lazy(xml_text, "metabolite")

            columns_to_concat = ["name", "synonyms", "iupac_names", "traditional_iupac"]
            # Apply the function row-wise to create a new column
            df['name'] = df.apply(concatCols, axis=1, colNames=columns_to_concat, sep='|')
            df = df.rename(_keys)
            if df is None:
                warnings.warn("Error during concatenation of columns in HMDB. Not Database created.")
                return 1
            writeData(df, db = "VMH")


if __name__ == '__main__':
    # Run as python -m src.dbhandling.reformatHMDB Databases/hmdb_metabolites.zip hmdb_metabolites.xml from the root folder
    # or as python  src/dbhandling/reformatHMDB.py Databases/hmdb_metabolites.zip hmdb_metabolites.xml
    main()

