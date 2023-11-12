def no_compartments(metabolites):
  # get metabolite identifiers without compartments for model metabolites
  # input: metabolites -> model.metabolites for model loaded with COBRApy
  # output: identifiers -> list of model metabolite identifiers without compartments
  identifiers = []
  for metabolite in metabolites:
    if metabolite.id[-1] == "]": # if compartments are [c]
      identifiers.append(metabolite.id[:-3])
    elif metabolite.id[-1] == ")": # if compartments are (c)
      identifiers.append(metabolite.id[:-3])
    elif metabolite.id[-2] == "_": # if compartments are _c
      identifiers.append(metabolite.id[:-2])
    elif metabolite.id[-3] == "_": # if compartments are _c0
      identifiers.append(metabolite.id[:-3])
  return identifiers