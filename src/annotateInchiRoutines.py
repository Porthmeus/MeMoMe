# Porthmeus
# 16.06.23
from typing import Optional

from src.matchMets import matchMetsByInchi, NeutraliseCharges
from rdkit import Chem, RDLogger
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
import warnings



# Define your function
def inchiToMol(inchi:str)->Chem.rdchem.Mol|None:
    if inchi is not None:
        return Chem.MolFromInchi(inchi)
    else:
        return None


# Define your function
def molToRDK(mol:Chem.rdchem.Mol)->ExplicitBitVect|None:
    if mol is not None:
        return Chem.RDKFingerprint(mol)
    else:
        return None


def molToNormalizedInchi(mol:Chem.rdchem.Mol, verbose = False) -> str|None:
    if verbose == False:
        RDLogger.DisableLog("rdApp.*")

    if mol is not None:
        return Chem.MolToInchi(mol)
    else:
        return None

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

def findOptimalInchi(inchis_: list[str], charge:int|None = None, verbose:bool = False) -> Optional[str]:
    """
    Find the "best" inchi string of equivalently annotated inchi strings.

    Sometimes it happens that the same molecule gets annotated with more than one inchi string. We need to handle this in some way. Usually it is just some form of isomer, which is not a big deal to take care of, but in case the Inchis are fundamentally different here we set following rules.
    1. Check which inchis match the charge of the molecule
    2. Use the inchi which has most trues after src.matchInchi.matchInchi() (this is the one which is most congruent with any other inchis in the list)
    3. If there are equals - use the most complex InChI string, that is the one with the most "/" in the string
    4. If this number is also equal, simply use the longest inchi string
    5. If this is also equal, use the first entry by chance

    Variables
    ---------
    inchis - list of inchi strings
    charge - integer denoting the charge of the metabolite as given by the model (!)
    """
    

    # first validate the inchi strings
    inchis = [x for x in inchis_ if validateInchi(x, verbose = verbose)]
    # return the inchi if there is only one left
    if len(inchis) == 1:
        return(inchis[0])

    # create a list of rdkit mol from inchis
    mols = [inchiToMol(x) for x in inchis]
    
    # rule1: check the charges
    if charge != None:
        inchis_correct_charge = [x for x,y in zip(inchis,mols) if Chem.GetFormalCharge(y) == charge]
        # if none of the inchis have the correct charge, keep them all and proceed
        if len(inchis_correct_charge) > 0:
            inchis = inchis_correct_charge
    
    # return the inchi if there is only one left
    if len(inchis) == 1:
        return(inchis[0])

    # rule 2
    matches = []
    for i in range(len(inchis)):
        k = 0
        for j in range(len(inchis)):
            if j != i:
                m1 = mols[i]
                m2 = mols[j]

                nminchi1 = molToNormalizedInchi(m1)
                nminchi2 = molToNormalizedInchi(m2)
                nt_m1 = NeutraliseCharges(m1)
                nt_m2 = NeutraliseCharges(m2)
                k = k + int(matchMetsByInchi(nminchi1, nminchi2, m1, m2, molToRDK(m1), molToRDK(m2), nt_m1, nt_m2, verbose = verbose)[0])
        matches.append(k)
    if len(matches) == 0:
        return None
    max_match = max(matches)
    inchis = [x for x, y in zip(inchis, matches) if y == max_match]

    # return the inchi if there is only one left
    if len(inchis) == 1:
        return(inchis[0])

    # rule 3
    counts = []
    for inchi in inchis:
        counts.append(inchi.count("/"))
    if len(counts) == 0:
        return None
    max_count = max(counts)
    inchis = [x for x, y in zip(inchis, counts) if y == max_count]

    # return the inchi if there is only one left
    if len(inchis) == 1:
        return(inchis[0])

    # rule 4
    counts = []
    for inchi in inchis:
        counts.append(len(inchi))
    if len(counts) == 0:
        return None
    max_count = max(counts)
    inchis = [x for x, y in zip(inchis, counts) if y == max_count]
    if len(inchis) == 0:
        return None

    # return the inchi if there is only one left
    if len(inchis) == 1:
        return(inchis[0])

    # rule 5. For determinism, we sort them first so we always get the same result
    inchis.sort()
    return inchis[0]



