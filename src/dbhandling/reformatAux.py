# Porthmeus
# 22.05.25
from src.download_db import get_config, get_database_path
import pandas as pd
import os
import json

def __json_to_dataframe(path) -> pd.DataFrame:
  with open(path, "r") as f: 
    vmh = json.load(f)["results"]
    return(pd.DataFrame(vmh))

def getData(db:str) -> pd.DataFrame:
    # loads database and the indentifier prefixes and returns them
    # get the database
    config = get_config()
    dbs_csv = ["BiGG","ModelSeed"]
    dbs_json = ["VMH"]
    if db in dbs_csv:
        dat = pd.read_csv(os.path.join(get_database_path(),config["databases"][db]["file"]),
                          sep = "\t",
                          low_memory=False)
    elif db in dbs_json:
        file = os.path.join(get_database_path(),config["databases"][db]["file"])
        dat = __json_to_dataframe(file)
    else:
        raise ValueError("db must be one of " + str(dbs_csv + dbs_json))
    return(dat)


def writeData(dat:pd.DataFrame, db:str) -> None:
    config = get_config()
    outfile = os.path.join(get_database_path(), config["databases"][db]["reformat"])
    dat.to_csv(outfile)
