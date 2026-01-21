from urllib.request import urlopen
import json
import warnings
from io import StringIO
import os
import pandas as pd
from src.download_db import get_config, get_database_path
from src.MeMoMetabolite import MeMoMetabolite


def origin_databases(metabolites:list[MeMoMetabolite]) -> dict:
    # identify origin databases for model metabolites
    # input: metabolites -> a list of MeMoMetabolites 
    # output: metabolite_namespace -> dictionary with fraction of model metabolites found in each database
    metabolite_namespace  = {}
    # convert the list of model metabolite identifiers without compartments to a dataframe
    df_metabolites = pd.DataFrame([met._id for met in metabolites],columns = ["model_metabolites"])
    # get the config file
    config = get_config()
    for db in config["databases"]:
      if db == "Identifiers":
        continue
      db_path =  os.path.join(get_database_path(), config["databases"][db]["reformat"])
      try:
          df = pd.read_table(db_path, sep = ",")
          metabolites_count = df_metabolites["model_metabolites"].isin(df["id"]).sum()
          metabolite_namespace[db] = metabolites_count / len(df_metabolites["model_metabolites"])
      except FileNotFoundError as e:
          warnings.warn(str(e))

    return metabolite_namespace
