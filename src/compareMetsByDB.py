# Porthmeus
# 02.11.23

from src.MeMoMetabolite import MeMoMetabolite
import pandas as pd

def compareMetsByDB(met1:MeMoMetabolite, met2:MeMoMetabolite) -> float:
    # use the database annotations to compare two metabolites
    # get common database entries
    commonDB = list(set(met1.annotations.keys()) & set(met2.annotations.keys()))
    # get the individual entries for each data base
    set1 = set([x+"."+item for x in commonDB for item in met1.annotations[x]])
    set2 = set([x+"."+item for x in commonDB for item in met2.annotations[x]])
    # calculate jaccard index
    inter_len = len(list(set1 & set2))
    union_len = len(list(set1 | set2))
    return(inter_len/union_len)
