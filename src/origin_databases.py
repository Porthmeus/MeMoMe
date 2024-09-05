from urllib.request import urlopen
import json
import warnings
from io import StringIO
import os
import pandas as pd
from src.download_db import get_config, get_database_path
from src.MeMoMetabolite import MeMoMetabolite


def VMH2DataFrame():
    # get the necessary metabolite info from VMH, and return it in a dataframe
    config = get_config()
    db_path =  os.path.join(get_database_path(), config["databases"]["VMH"]["file"])

    # load VMH as dataframe
    try:
        with open(db_path, "r") as json_file:
            VMH_json = json.load(json_file)
            VMH_json = json.dumps(VMH_json['results'])
            VMH = pd.read_json(StringIO(VMH_json))
    except FileNotFoundError as e:
        warnings.warn(e)
        VMH = pd.DataFrame()

    # return the file as dataframe
    return VMH



def origin_databases(metabolites:list[MeMoMetabolite]) -> dict:
    # identify origin databases for model metabolites
    # input: metabolites -> a list of MeMoMetabolites 
    # output: metabolite_namespace -> dictionary with fraction of model metabolites found in each database

    # convert the list of model metabolite identifiers without compartments to a dataframe
    df_metabolites = pd.DataFrame([met._id for met in metabolites],columns = ["model_metabolites"])

    # get the config file
    config = get_config()

    # get metabolite info from BiGG
    db_path =  os.path.join(get_database_path(), config["databases"]["BiGG"]["file"])
    try:
        df_BiGG = pd.read_table(db_path)
    except FileNotFoundError as e:
        warnings.warn(e)
        df_Bigg = pd.DataFrame()
    # find number of model metabolites in BiGG database
    metabolitesBiGG_count = df_metabolites["model_metabolites"].isin(df_BiGG["universal_bigg_id"]).sum()

    # get metabolite info from VMH
    df_VMH = VMH2DataFrame()
    # find number of model metabolites in VMH database
    metabolitesVMH_count = df_metabolites["model_metabolites"].isin(df_VMH["abbreviation"]).sum()

    # get metabolite info from ModelSEED
    db_path =  os.path.join(get_database_path(), config["databases"]["ModelSeed"]["file"])
    try:
        df_ModelSEED = pd.read_table(db_path, low_memory=False)
    except FileNotFoundError as e:
        warnings.warn(e)
        df_ModelSEED = pd.DataFrame()
    # find number of model metabolites in ModelSEED database
    metabolitesModelSEED_count = df_metabolites["model_metabolites"].isin(df_ModelSEED["id"]).sum()

    # dictionary with percent of model metabolites found in each database
    metabolite_namespace = {"BiGG": metabolitesBiGG_count / len(df_metabolites["model_metabolites"]),
                          "VMH": metabolitesVMH_count / len(df_metabolites["model_metabolites"]),
                          "ModelSEED": metabolitesModelSEED_count / len(df_metabolites["model_metabolites"])}
    return metabolite_namespace
