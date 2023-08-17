import pandas as pd



def ModelSEED():

  # load ModelSEED as dataframe
  ModelSEED = pd.read_csv("https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/compounds.tsv", sep="\t", low_memory=False).fillna("")

  # get final dataframe with info from ModelSEED
  df_ModelSEED = ModelSEED[["id","name","formula","mass","inchikey","charge","deltag","deltagerr","pka","pkb","aliases","smiles"]]
  return df_ModelSEED