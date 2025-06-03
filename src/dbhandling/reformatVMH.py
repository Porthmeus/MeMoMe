# Porthmeus
# 10.04.25

# libraries
from src.dbhandling.reformatAux import *
from src.annotation.annotateInchiRoutines import *
import pandas as pd


def getAnnosPerEntry(dat:pd.DataFrame, met:str) -> dict[str,list[str]]:
    """ takes the dataframe for vmh entries and a str from the metabolite id and returns an annotation dictionary for selected annotations in VMH """
    keys = {"biggId":"bigg.metabolite",
            "keggId":"kegg.compound",
            "chemspider":"chemspider",
            "cheBlId":"chebi",
            "biocyc":"biocyc", # not supported by identifiers.org
            "food_db":"foodb.compound",
            "hmdb":"hmdb",
            "iuphar_id":"iuphar.ligand",
            "metanetx":"metanetx", # not supported by identifiers.org
            "metlin":"metlin",
            "fda_id":"unii",
            "pubChemId":"pubchem.compound",
            "seed":"seed.compound"}
     
    dat_sel = dat.loc[dat.abbreviation == met,]
    anno = {"vmhmetabolite":[met]}
    for key in keys.keys():
        val = dat_sel[key].iloc[0]
        if val != None:
            anno.setdefault(keys[key], []).append(val)
    return(anno)

def getAnnos(dat:pd.DataFrame) -> pd.Series:
    """takes the data frame from VMH request and returns a pandas series of strings. Each string is a dictionary for the database annotations and can be simply evaluated (eval()) to be transformed in the correct format"""
    anno_series = pd.Series()
    for i in range(len(dat)):
        met = dat.iloc[i]["abbreviation"]
        anno = getAnnosPerEntry(dat,met)
        anno_series[len(anno_series)] = str(anno)
    anno_series.index = list(dat["abbreviation"])
    return(anno_series)

def concatNames(x:pd.Series) -> str:
    """ handles the concatenation of the different fields which we consider a name of a metabolite. Will return a single string containing the names which are seperated by a '|'."""
    iupac = x["iupac"]
    name = x["fullName"]
    alias = x["synonyms"]
    concat = []
    if name != "" and name != None:
        concat.append(name)
    if iupac != "" and iupac != None:
        concat.append(iupac)
    if alias != "" and alias != None:
        alias = alias.replace("***","|")
        concat.append(alias)
    concat = "|".join(concat)
    return(concat)
# for names check fullName iupac and alias
# alias contains a list of names seperated by "***" 

def reformatVMH() -> None:
    """ Takes the VMH database, reformats it and writes a csv with standard columns (id, inchi, name, DBs) + additional information to disk """

    vmh = getData("VMH")
    vmh.index = vmh.abbreviation
    # get names, inchis and database annotations
    names = vmh.apply(concatNames, axis = 1).fillna(value = "")
    inchis = vmh.inchiString.fillna(value = "")
    dbs = getAnnos(vmh)

    # concat the different series/dataframes
    dat_all = pd.concat([names,inchis,dbs], axis = 1)
    dat_all.columns = ["name","inchi","DBs"]
    dat_all["id"] = dat_all.index
    dat_all = dat_all.loc[:,["id","name","inchi","DBs"]]
    dat_all = pd.concat([dat_all, vmh], axis = 1)
    
    # save data
    writeData(dat_all, db = "VMH")
    
