import pandas as pd
from ast import literal_eval



def BiGG():
  # get the necessary metabolite info from BiGG, and return it in a dataframe

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

  # load the file with InChIKey, formulae and charges for metabolites from BiGG API as a dataframe
  data_df = pd.read_csv("BiGG_API.csv").fillna("")
  # convert dataframe entries to lists
  data_df["formulae"] = data_df["formulae"].apply(literal_eval)
  data_df["charges"] = data_df["charges"].apply(literal_eval)
  # keep formulae and charges from BiGG API
  BiGG["formula"] = data_df["formulae"] # each row is a list
  BiGG["charge"] = data_df["charges"] # each row is a list

  # keep existing InChIKey for metabolites
  BiGG["inchikey"] = BiGG["database_links"].str.split("InChI Key: https://identifiers.org/inchikey/").str[1].str.split("; ").str[0].fillna("")

  # get final dataframe with info from BiGG
  df_BiGG = BiGG[["universal_bigg_id","name","formula","charge","inchikey"]]
  return df_BiGG