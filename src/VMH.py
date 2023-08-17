from urllib.request import urlopen
import json
from io import StringIO



def VMH():

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
  return df_VMH