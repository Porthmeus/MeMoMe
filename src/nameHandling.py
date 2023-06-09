# Porthmeus
# 09.06.23

# function to handle some specifics of the naming convention in SBML files and metabolic modeling


def removeMetabolitePrefixSuffix(istr: str) -> str:
    """
    Takes a string and removes the prefix (M_) and the compartment string.
    Returns a cleaned string
    """
    ostr = re.sub("M_(.*)_[a-z]\d?$","\g<1>",istr)
    return(ostr)

