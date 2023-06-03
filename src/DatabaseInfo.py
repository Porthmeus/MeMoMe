import pandas as pd
from ast import literal_eval
from urllib.request import urlopen
import json
from io import StringIO



# BiGG

# load BiGG as dataframe
BiGG = pd.read_csv("http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt", sep="\t").fillna("")

# load the file with InChIKey, formulae and charges for metabolites from BiGG API as a dataframe
data_df = pd.read_csv("BiGG_API.csv").fillna("")

# convert dataframe entries to lists
data_df["InChIKey"] = data_df["InChIKey"].apply(lambda x: literal_eval(x) if x != "" else x)
data_df["formulae"] = data_df["formulae"].apply(literal_eval)
data_df["charges"] = data_df["charges"].apply(literal_eval)

# skip names that don't exist as metabolites in BiGG API
BiGG = BiGG[~BiGG["universal_bigg_id"].isin(["Recon3D", "iYS854"])].reset_index()

# keep InChIKey from BiGG API
BiGG["inchikey"] = data_df["InChIKey"].apply(lambda x: x[0] if x != "" else x)

# keep formulae and charges from BiGG API
BiGG["formula"] = data_df["formulae"] # each row is a list
BiGG["charge"] = data_df["charges"] # each row is a list

# get final dataframe with info from BiGG
df_BiGG = BiGG[["universal_bigg_id","name","formula","charge","inchikey"]]



# VMH

# load VMH as dataframe

url = urlopen("https://www.vmh.life/_api/metabolites/?format=jsonp&page_size=9999999&_dc=1670601370202&page=1&start=0&limit=9999999&callback=Ext.data.JsonP.callback19")

json_file = url.read().decode("utf-8")

json_file = json_file[len("Ext.data.JsonP.callback19("):]
json_file = json_file[:-len(");")]

VMH_json = json.loads(json_file)

VMH_json = json.dumps(VMH_json['results'])

VMH = pd.read_json(StringIO(VMH_json))

# get final dataframe with info from VMH
df_VMH = VMH[["abbreviation","fullName","synonyms","neutralFormula","chargedFormula","charge","avgmolweight","monoisotopicweight","inchiString","inchiKey","smile"]]



# ModelSEED

# load ModelSEED as dataframe
ModelSEED = pd.read_csv("https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/compounds.tsv", sep="\t", low_memory=False).fillna("")

# get final dataframe with info from ModelSEED
df_ModelSEED = ModelSEED[["id","name","formula","mass","inchikey","charge","deltag","deltagerr","pka","pkb","aliases","smiles"]]
