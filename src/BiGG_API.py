import pandas as pd
import requests
import time

# load BiGG as dataframe
BiGG = pd.read_csv("http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt", sep="\t").fillna("")

# get InChIKey, formulae and charges for metabolites from BiGG API
data_list = []
for metabolite in BiGG["universal_bigg_id"]:
  endpoint = f"http://bigg.ucsd.edu/api/v2/universal/metabolites/{metabolite}"
  r = requests.get(endpoint)
  if not r.ok: # skip names that don't exist as metabolites  ### BiGG = BiGG[~BiGG["universal_bigg_id"].isin(["Recon3D", "iYS854"])]
    continue
  result = r.json()
  if "InChI Key" in result["database_links"]:
    InChIKey = [inchikey["id"] for inchikey in result["database_links"]["InChI Key"]]
  else:
    InChIKey = ""
  formulae = result["formulae"]
  charges = result["charges"]
  row = [InChIKey,  formulae, charges] # each one is a list with strings
  data_list.append(row)
  time.sleep(0.1) # delay some seconds

# save the list with the retrived info as a dataframe and in a file
data_df = pd.DataFrame(data_list, columns = ["InChIKey","formulae","charges"])
data_df.to_csv("BiGG_API.csv",index=False)