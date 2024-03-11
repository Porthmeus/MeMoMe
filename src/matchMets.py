# Porthmeus
# 02.11.23

#from src.MeMoMetabolite import MeMoMetabolite
from __future__ import annotations
import pandas as pd
from src.matchInchiRoutines import *
from rdkit import Chem, RDLogger
from Levenshtein import ratio
from warnings import warn

def matchMetsByInchi(inchi1:str, inchi2:str, mol1: Chem.rdchem.Mol, mol2: Chem.rdchem.Mol, fp1, fp2, verbose:bool = False) -> tuple:
    ''' A combination of different matching algorithms for the Inchi-Strings, which end in a simple yes/no answer, whether the molecules are the same'''

    # turn off chattiness of rdkit
    if verbose == False:
        RDLogger.DisableLog("rdApp.*")

    # rdkit is somewhat chatty. To supress the warning:
    # check whether the molecules are differently charged, if so neutralize the charges
#try:
    m1 = mol1
    m2 = mol2
    charge1 = Chem.GetFormalCharge(m1)
    charge2 = Chem.GetFormalCharge(m2)
    chargeDiff = charge1 - charge2
    if charge1 != charge2:
        m1 = NeutraliseCharges(m1)
        m2 = NeutraliseCharges(m2)
        inchi1 = Chem.MolToInchi(m1)
        inchi2 = Chem.MolToInchi(m2)
# except:
#     warn("Something went wrong with the InChI comparison\nmet1: " + inchi1 + "\nmet2: " + inchi2)
#     same = False
#     chargeDiff = None
# else:
    # check the fingerprints - they should be the same, if not, check at the same time, whether after H removal there is only one atom left.
    same = False
    if compareInchiByFingerprint0(fp1, fp2, verbose = verbose) == 1.0 or m1.GetNumAtoms() == m2.GetNumAtoms() == 1:

        # checking the canonical InchiString
        same = compareInchiByString0(m1, m2, verbose = verbose)
        if same == False:

            # if that is false check if there is one of the stereoisomers the
            # same
            same = compareInchiByStereoIsomer0(m1,m2, verbose = verbose)
# finally:
    return(same, chargeDiff)

def matchMetsByDB(met1:MeMoMetabolite, met2:MeMoMetabolite) -> float:
    # use the database annotations to compare two metabolites
    # get common database entries
    commonDB = list(set(met1.annotations.keys()) & set(met2.annotations.keys()))
    # get the individual entries for each data base
    set1 = set([x+"."+item for x in commonDB for item in met1.annotations[x]])
    set2 = set([x+"."+item for x in commonDB for item in met2.annotations[x]])
    # calculate jaccard index
    inter_len = len(list(set1 & set2))
    union_len = len(list(set1 | set2))
    if union_len == 0:
        sim = 0
    else: 
        sim = inter_len/union_len
    return(sim)

def matchMetsByName(met1:MeMoMetabolite, met2:MeMoMetabolite) -> float:
    # use the metabolite names to match
    score = 0
    for name1 in met1.names:
        for name2 in met2.names:
            tmp_score = ratio(name1,name2)
            if tmp_score > score:
                score = tmp_score
    return(score)
    
    

