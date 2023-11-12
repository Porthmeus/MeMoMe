import pandas as pd
import requests
import time
from tqdm import tqdm
import os
from pathlib import Path


def BiGG_API(output = Path("BiGG_API.csv")) -> None:
  # extract metabolite formulae and charges from BiGG API, and save them in BiGG_API.csv

  # load BiGG as dataframe
  BiGG = pd.read_csv("http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt", sep="\t").fillna("")

  # fix the 3 wrongly parsed metabolites in BiGG
  # find their indeces
  i1 = BiGG[BiGG["universal_bigg_id"] == "Recon3D"].index
  i2 = BiGG[BiGG["universal_bigg_id"] == "iYS854"].index
  i = i1.append(i2)
  # shuffle their information at the right columns
  temp = BiGG["universal_bigg_id"][i]
  BiGG["name"][i] = BiGG["bigg_id"][i]
  BiGG["bigg_id"][i] = ""
  BiGG["universal_bigg_id"][i] = BiGG["model_list"][i].str.split("; ").str[0]
  BiGG["old_bigg_ids"][i] = BiGG["model_list"][i].str.split("; ").str[1]
  BiGG["model_list"][i] = temp

  if os.path.exists(output):
      bigg_api = pd.read_csv(output)
      already_in = list(bigg_api.metabolite)
  else:
      already_in = []

  # get formulae and charges for metabolites from BiGG API
  data_list = []
  for metabolite in tqdm(BiGG["universal_bigg_id"].unique()):
    if metabolite not in already_in:
        endpoint = f"http://bigg.ucsd.edu/api/v2/universal/metabolites/{metabolite}"
        r = requests.get(endpoint)
        result = r.json()
        formulae = result["formulae"]
        charges = result["charges"]
        df = pd.DataFrame({"metabolite":[metabolite],
            "formulae": [formulae],
            "charges" : [charges]}) # each one is a list with strings
        if os.path.exists(output):
            with open(output, "a") as ff:
                df.to_csv(ff, index = False)
        else:
            df.to_csv(output, index = False)
        #data_list.append(row)
        time.sleep(0.1) # delay some seconds

