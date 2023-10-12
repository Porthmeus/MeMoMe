#SteFlor, adapted from Porthmeus and DafniGi
# 12/10/2023

import re
from typing import Optional
import sys
sys.path.append('../')
from src.MeMoMetabolite import MeMoMetabolite

def handle_metabolites_prefix_suffix(metabolites: list[MeMoMetabolite]): -> list:
    '''
    Parses the SBML model object, and extracts metabolite after prefix and suffix removal 
    Parameters
    ----------
        sbml_file_url : pathlib.Path -  an sbml model file path
    '''
    identifiers = []
    for met in metabolites:
        string = met.id
        if "cpd" in string: # extract cpd
	    cpd_index = string.find("cpd")
	    met_id = "cpd"
	    # go to the index of 'cpd' + 3(the length) of cpd because we now expect
	    # the numerical id
	    for i in range(cpd_index + 3, len(string)):
	        # only append the digits and stop as soon as you encounter a non-digit
	        if string[i].isdigit():
	            met_id += string[i]
	        else:
	            break
        else: # remove 'M_' prefix and compartment string from metabolite id, preserving content between them.
	    met_id = re.sub("^M_(.*)_[a-z]\d?$","\g<1>", string)

        identifiers = identifiers.append(met_id)
    return identifiers
