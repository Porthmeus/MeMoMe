# Porthmeus
# 02.11.23

#from src.MeMoMetabolite import MeMoMetabolite
from __future__ import annotations
import pandas as pd
from src.annotateInchiRoutines import *
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from Levenshtein import ratio
from warnings import warn
from collections import namedtuple

def matchMetsByInchi(nminchi1: str,
        nminchi2: str,
        mol1: Chem.rdchem.Mol,
        mol2: Chem.rdchem.Mol,
        charge_neutralized_inchi1:str,
        charge_neutralized_inchi2:str,
        charge1:int,
        charge2:int,
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
    # check for charge differences and correct them if necessary
    chargeDiff = charge1 - charge2
    if charge1 != charge2:
        nminchi1 = charge_neutralized_inchi1
        nminchi2 = charge_neutralized_inchi2
    
    # check if we already have the same inchi string - normalized inchis also rectify tautorism
    same = nminchi1 == nminchi2
    

    # if that is not the case check if we have the same molecule but different stereoisomers
    if same == False:
        # remove the layer for cis-trans isomerism and check again
        nminchi1 = "/".join([x for x in nminchi1.split("/") if not x.startswith("b")])
        nminchi2 = "/".join([x for x in nminchi2.split("/") if not x.startswith("b")])
        same = nminchi1 == nminchi2
        #print(nminchi1, nminchi2)
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
    
    

