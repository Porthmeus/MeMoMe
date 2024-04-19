# Porthmeus
# 02.11.23

#from src.MeMoMetabolite import MeMoMetabolite
from __future__ import annotations
import pandas as pd
from src.matchInchiRoutines import *
from rdkit import Chem, RDLogger
from Levenshtein import ratio
from warnings import warn
from collections import namedtuple

def matchMetsByInchi(nminchi1: str,
        nminchi2: str,
        mol1: Chem.rdchem.Mol,
        mol2: Chem.rdchem.Mol,
        fp1,
        fp2,
        charge_neutralized_mol1,
        charge_neutralized_mol2,
        verbose:bool = False) -> tuple:
    ''' A combination of different matching algorithms for the Inchi-Strings, which end in a simple yes/no answer, whether the molecules are the same'''

    # rdkit is somewhat chatty. To supress the warning:
    # turn off chattiness of rdkit
    if verbose == False:
        RDLogger.DisableLog("rdApp.*")

    # check whether the molecules are differently charged, if so neutralize the charges
#try:
    m1 = mol1
    m2 = mol2
    charge1 = Chem.GetFormalCharge(m1)
    charge2 = Chem.GetFormalCharge(m2)
    chargeDiff = charge1 - charge2
    if charge1 != charge2:
        m1 = charge_neutralized_mol1
        m2 = charge_neutralized_mol2
 #       m1 = NeutraliseCharges(m1)
 #       m2 = NeutraliseCharges(m2)
    
    same = False
    # check the fingerprint
    if compareInchiByFingerprint0(fp1, fp2, verbose = verbose) == 1.0 or m1.GetNumAtoms() == m2.GetNumAtoms() == 1:
        same = nminchi1 == nminchi2
        if same == False:
            # if that is false check if there is one of the stereoisomers the
            # same
            same = compareInchiByStereoIsomer0(m1,m2, verbose = verbose)
# finally:
    return(same, chargeDiff)



DBResult = namedtuple("DBResult", ["commonDBs", "commonIds", "allIds", "score"])

def matchMetsByDB(met1:MeMoMetabolite, met2:MeMoMetabolite) -> DBResult:
    # use the database annotations to compare two metabolites
    # get common database entries
    commonDB = list(set(met1.annotations.keys()) & set(met2.annotations.keys()))
    # get the individual entries for each data base
    set1 = set([x+"."+item for x in commonDB for item in met1.annotations[x]])
    set2 = set([x+"."+item for x in commonDB for item in met2.annotations[x]])
    # calculate jaccard index
    inter = list(set1 & set2)
    unio = list(set1 | set2)
    inter_len = len(inter)
    union_len = len(unio)
    if union_len == 0:
        sim = 0
    else: 
        sim = inter_len/union_len
    # TODO: This sort is only so the result is always the samr order which is
    # necessary for the tests. If this slows down the program, think of another soluion
    inter.sort()
    unio.sort()
    return(DBResult(commonDB, str(inter), str(unio), sim))


NamedResult = namedtuple("NamedResult", ["name_id1", "name_id2", "score"])

def matchMetsByName(met1:MeMoMetabolite, met2:MeMoMetabolite) -> NamedResult:
    # use the metabolite names to match
    score = 0
    # Name with the highest score in model 1
    temp_name_1 = ""
    # Name with the highest score in model 2:w
    temp_name_2 = ""
    for name1 in met1.names:
        for name2 in met2.names:
            tmp_score = ratio(name1,name2)
            if tmp_score > score:
                temp_name_1 = name1
                temp_name_2 = name2
                score = tmp_score
    return (NamedResult(temp_name_1, temp_name_2, score))
    
    

