# Porthmeus
# 27.04.23

from src.matchInchiRoutines import *
from rdkit import Chem

def matchInchi(inchi1:str, inchi2:str) -> tuple:
    ''' A combination of different matching algorithms for the Inchi-Strings, which end in a simple yes/no answer, whether the molecules are the same'''

    # check whether the molecules are differently charged, if so neutralize the charges
    m1 = Chem.MolFromInchi(inchi1)
    m2 = Chem.MolFromInchi(inchi2)
    charge1 = Chem.GetFormalCharge(m1)
    charge2 = Chem.GetFormalCharge(m2)
    chargeDiff = charge1 - charge2
    if charge1 != charge2:
        m1 = NeutraliseCharges(m1)
        m2 = NeutraliseCharges(m2)
        inchi1 = Chem.MolToInchi(m1)
        inchi2 = Chem.MolToInchi(m2)

    # check the fingerprints - they should be the same, if not, check at the same time, whether after H removal there is only one atom left.
    same = False
    if compareInchiByFingerprint(inchi1,inchi2) == 1.0 or m1.GetNumAtoms() == m2.GetNumAtoms() == 1:

        # checking the canonical InchiString
        same = compareInchiByString(inchi1, inchi2)
        if same == False:

            # if that is false check if there is one of the stereoisomers the
            # same
            same = compareInchiByStereoIsomer(inchi1,inchi2)
    return(same, chargeDiff)
