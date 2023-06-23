# Porthmeus
# 09.06.23

# function to handle some specifics of the naming convention in SBML files and metabolic modeling

import re
from typing import Optional

def removeMetabolitePrefixSuffix(istr: str) -> str:
    """
    Takes a string and removes the prefix (M_) and the compartment string.
    Returns a cleaned string
    """
    ostr = re.sub("^M_(.*)_[a-z]\d?$","\g<1>",istr)
    return(ostr)


def extract_cpd(string: str) -> Optional[str]:
    """
    :param string: The string to extract the cpd from
    :return: Returns the extract cpd or none if the given string did not contain `cpd` as a substring
    """
    cpd_index = string.find("cpd")
    # check if the string contrains 'cpd'
    if(cpd_index == -1):
        return None
    cpd_id = "cpd"
    # if it does, go to the index of 'cpd' + 3(the length) of cpd because we now expect
    # the numerical id
    for i in range(cpd_index + 3, len(string)):
        # only append the digits and stop as soon as you encounter a non-digit
        if string[i].isdigit():
            cpd_id += string[i]
        else:
            break
    return cpd_id
