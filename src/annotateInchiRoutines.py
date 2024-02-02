# Porthmeus
# 16.06.23

from src.matchMets import matchMetsByInchi
from rdkit import Chem, RDLogger
import warnings

def validateInchi(inchi:str, verbose:bool = False) -> bool:
    """
    Takes an inchi string and checks whether an inchi is correct or whether it contains mistakes which are not handled by rdkit
    """
    # silence rdkit
    if verbose == False:
        RDLogger.DisableLog("rdApp.*")

    # standardize inchi and at the same time validate the string
    try:
        m = Chem.MolFromInchi(inchi)
        log = Chem.SanitizeMol(m)
        inchi = Chem.MolToInchi(m)
        validated = True
    except:
        validated = False
    return(validated)

def smile2inchi(smile:str, verbose:bool = False) -> str:
    """ Takes a smile and returns inchi_string """
    
    # silence rdkit
    if verbose == False:
        RDLogger.DisableLog("rdApp.*")

    # load smiles and convert to inchi
        m = Chem.MolFromSmiles(smile)
    try:
        log = Chem.SanitizeMol(m)
    except:
        warnings.warn("Could not sanitize {mol}".format(mol = smile))
    inchi = Chem.MolToInchi(m)
    return(inchi)

def findOptimalInchi(inchis: list[str], verbose:bool = False) -> str:
    """
    Find the "best" inchi string of equivalently annotated inchi strings.

    Sometimes it happens that the same molecule gets annotated with more than one inchi string. We need to handle this in some way. Usually it is just some form of isomer, which is not a big deal to take care of, but in case the Inchis are fundamentally different here we set following rules.
    1. Use the inchi which has most trues after src.matchInchi.matchInchi()
    2. If there are equals - use the most complex InChI string, that is the one with the most "/" in the string
    3. If this number is also equal, simply use the longest inchi string
    4. If this is also equal, use the first entry by chance
    """
    

    # first validate the inchi strings
    inchis = [x for x in inchis if validateInchi(x, verbose = verbose)]
    # return the inchi if there is only one left
    if len(inchis) == 1:
        return(inchis[0])

    # rule 1
    matches = []
    for i in range(len(inchis)):
        k = 0
        for j in range(len(inchis)):
            if j != i:
                k = k + int(matchMetsByInchi(inchis[i], inchis[j], verbose = verbose)[0])
        matches.append(k)

    max_match = max(matches)
    inchis = [x for x, y in zip(inchis, matches) if y == max_match]

    # rule 2
    counts = []
    for inchi in inchis:
        counts.append(inchi.count("/"))
    max_count = max(counts)
    inchis = [x for x, y in zip(inchis, counts) if y == max_count]

    # rule 3
    counts = []
    for inchi in inchis:
        counts.append(len(inchi))
    max_count = max(counts)
    inchis = [x for x, y in zip(inchis, counts) if y == max_count]

    # rule 4. For determinism, we sort them first so we always get the same result
    inchis.sort()
    return inchis[0]
