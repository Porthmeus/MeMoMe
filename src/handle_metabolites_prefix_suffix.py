#SteFlor, adapted from Porthmeus and DafniGi
# 12/10/2023

import re
from typing import Optional
import sys
sys.path.append('../')

def handle_metabolites_prefix_suffix(met_id: str) -> str:
    '''
    # function to handle some specifics of the naming convention in SBML files and metabolic modeling
    ----------
        met_id : str -  the id to be manipulated
        identifier: str - the id without prefix and suffix
    '''

    identifier = met_id
    if "cpd" in met_id: # extract cpd
        cpd_index = met_id.find("cpd")
        identifier = "cpd"
        # go to the index of 'cpd' + 3(the length) of cpd because we now expect
        # the numerical id
        for i in range(cpd_index + 3, len(met_id)):
            # only append the digits and stop as soon as you encounter a non-digit
            if met_id[i].isdigit():
                identifier += met_id[i]
            else:
                break
    else: # remove 'M_' prefix and compartment string from metabolite id, preserving content between them.
        identifier = re.sub("^M_(.*)_[a-z]\d?$","\g<1>", met_id)
        identifier = re.sub("^M_(.*)__91__[a-z]\d?__93__$","\g<1>", identifier)
        identifier = re.sub("^M_(.*)__40__[a-z]\d?__41__$","\g<1>", identifier)
    return identifier
