# SteFlor, adapted from Porthmeus and DafniGi
# 12/10/2023

import re
from typing import Optional
import sys

sys.path.append("../")


def handle_metabolites_prefix_suffix(met_id: str) -> str:
    """
    # function to handle some specifics of the naming convention in SBML files and metabolic modeling
    ----------
        met_id : str -  the id to be manipulated
        identifier: str - the id without prefix and suffix
    """

    identifier = met_id
    if "cpd" in met_id.lower():
        if (
            "cpd" not in met_id
        ):  # the "cpd" part of the string contains one or more uppercase charachters. It will be considered as an invalid id
            identifier = None
        else:  # extract cpd
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
    else:  # remove 'M_' prefix and compartment string from metabolite id, preserving content between them.
        identifier = re.sub(r"^(M_)?(.*)_[a-z]\d?$", r"\g<2>", met_id)
        identifier = re.sub(r"^(M_)?(.*)__91__[a-z]\d?__93__$", r"\g<2>", identifier)
        identifier = re.sub(r"^(M_)?(.*)__40__[a-z]\d?__41__$", r"\g<2>", identifier)
        identifier = re.sub(r"^(M_)?(.*)\[[a-z]\d?\]$", r"\g<2>", identifier)
        identifier = re.sub(r"^(M_)?(.*)\([a-z]\d?\)$", r"\g<2>", identifier)
    return identifier


def validateInchi(inchi_string: str) -> str:
    """Checks if an inchi string conforms to some standards, if not, returns None"""
    correct_inchi = inchi_string
    if inchi_string is None or inchi_string == "":
        correct_inchi = None
    elif not inchi_string.startswith("InChI="):
        correct_inchi = None
    elif not "/" in inchi_string:
        correct_inchi = None
    elif len(inchi_string.split("/")) < 3:
        correct_inchi = None
    return correct_inchi
