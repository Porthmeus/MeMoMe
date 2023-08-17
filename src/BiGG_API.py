import pandas as pd
import requests
import time
from tqdm import tqdm



def BiGG_API():

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


  # get formulae and charges for metabolites from BiGG API
  data_list = []
  for metabolite in tqdm(BiGG["universal_bigg_id"]):
    endpoint = f"http://bigg.ucsd.edu/api/v2/universal/metabolites/{metabolite}"
    r = requests.get(endpoint)
    result = r.json()
    formulae = result["formulae"]
    charges = result["charges"]
    row = [formulae, charges] # each one is a list with strings
    data_list.append(row)
    time.sleep(0.1) # delay

  # save the list with the retrieved info as a dataframe and in a file
  data_df = pd.DataFrame(data_list, columns = ["formulae","charges"])
  data_df.to_csv("BiGG_API.csv",index=False)