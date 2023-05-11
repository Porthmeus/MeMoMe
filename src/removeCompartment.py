#Hadjar
#defined function and test cases
#10.05.2023

#Stefano
#11.05.2023
#added a few name suffixes, splitted into code and test files, splitted into one function for ids and one for names


import re

    """
    Removes the suffixes (which refer to the compartment in which the metabolites are found) from the metabolite ids, as well as the prefixes (which refer to the organism in which the metabolites are found) from the metabolite ids.

    Args:
        metabolite_ids (list): A list of metabolite ids from which the compartment suffixes, and the organism prefixes, need to be removed if present.

    Returns:
        list: A list of the metabolites ids after removal of the compartment suffixes and organism prefixes.
    """

    #TODO:  add handling of _comp0 suffixes

def removeCompAndOrganismFromMetIds(metabolite_ids: list):
    # Create a new list to store the metabolite ids without compartments
    new_metabolite_ids = []
    # Define a regular expression pattern to match the format "met['c']" or "met[c]"
    pattern = r"\w+\[?\w*\]?"
    # Loop through the original list of metabolite ids and remove the compartments
    for met_id in new_metabolite_ids:
        if met_id.startswith('M_'):
            new_met_id = met_id.split('_')[1]
            new_metabolite_ids.append(new_met_id)
        else:
            if re.fullmatch(pattern, met_id):
                new_met_id = met_id.split("[")[0]
                new_metabolite_ids.append(new_met_id)
            else:
                new_metabolite_ids.append(met_id)
    return new_metabolite_ids


#def removeCompAndOrganismFromMetNames(metabolite_names: list):

#    """
#    Removes the suffixes (which refer to the compartment in which the metabolites are found) from the metabolite names, as well as the prefixes (which refer to the organism in which the metabolites are found) from the metabolite names.

#    Args:
#        metabolite_names (list): A list of metabolite names from which the compartment suffixes, and the organism prefixes, need to be removed if present.

#    Returns:
#        list: A list of the metabolites names after removal of the compartment suffixes and organism prefixes.
#    """
#
#    # Create a new list to store the metabolite names without compartments
#    new_metabolite_names = []
#    # Define a regular expression pattern to match the format "met-c0"
#    pattern = r"\w+\-?\w*?"
#    # Loop through the original list of metabolite names and remove the compartments
#    for name in metabolite_names:
#        if re.fullmatch(pattern, name):
#            new_name = name.split("-")[0]
#            new_metabolite_names.append(new_name)
#        else:
#            new_metabolite_names.append(name)
#    return new_metabolite_names
