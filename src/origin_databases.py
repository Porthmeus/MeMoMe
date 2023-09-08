import pandas as pd



def origin_databases(metabolites):
  # identify origin databases for model metabolites
  # input: metabolites -> dataframe of model metabolite identifiers without compartments
  # output: metabolite_namespace -> dictionary with percent of model metabolites found in each database

  # convert the list of model metabolite identifiers without compartments to a dataframe
  df_metabolites = pd.DataFrame(metabolites,columns = ["model_metabolites"])

  # get metabolite info from BiGG
  df_BiGG = BiGG()
  # find number of model metabolites in BiGG database
  metabolitesBiGG_count = df_metabolites["model_metabolites"].isin(df_BiGG["universal_bigg_id"]).sum()

  # get metabolite info from VMH
  df_VMH = VMH()
  # find number of model metabolites in VMH database
  metabolitesVMH_count = df_metabolites["model_metabolites"].isin(df_VMH["abbreviation"]).sum()

  # get metabolite info from ModelSEED
  df_ModelSEED = ModelSEED()
  # find number of model metabolites in ModelSEED database
  metabolitesModelSEED_count = df_metabolites["model_metabolites"].isin(df_ModelSEED["id"]).sum()

  # dictionary with percent of model metabolites found in each database
  metabolite_namespace = {"BiGG": f'{(metabolitesBiGG_count / len(df_metabolites["model_metabolites"])) * 100:.1f}%',
                        "VMH": f'{(metabolitesVMH_count / len(df_metabolites["model_metabolites"])) * 100:.1f}%',
                        "ModelSEED": f'{(metabolitesModelSEED_count / len(df_metabolites["model_metabolites"])) * 100:.1f}%'}
  return metabolite_namespace