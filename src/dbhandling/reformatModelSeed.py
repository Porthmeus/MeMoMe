# Porthmeus
# 17.01.25

from src.download_db import get_config, get_database_path
from src.annotation.annotateInchiRoutines import smile2inchi
from io import StringIO
import pandas as pd
import os
import json

# this needs to hang in global space, otherwise we need to read the info over and over again
try: 
    # load the prefixes for the identifies.org list that is needed for the
    # correction of the key in the annotation dictionary
    config = get_config()
    identifiers_file = os.path.join(get_database_path(), config["databases"]["Identifiers"]["file"])
    with open(identifiers_file, "r") as json_file: 
        identifiers = json.load(json_file)
except FileNotFoundError as e:
    warnings.warn(str(e))
    identifiers = None
    # Rethrow exception because we want don't allow missing dbs
    # not working here, but there will be check further down again
if identifiers != None:
    identifier_prefixes = json.dumps(identifiers["_embedded"]["namespaces"])
    identifier_prefixes = pd.read_json(StringIO(identifier_prefixes))
    identifier_prefixes = list(identifier_prefixes["prefix"])

def getData() -> pd.DataFrame:
    # loads the modelSeed database and the indentifier prefixes and returns them
    # get the database
    config = get_config()
    dat = pd.read_csv(os.path.join(get_database_path(),config["databases"]["ModelSeed"]["file"]),
                      sep = "\t",
                      low_memory=False)
    return(dat)

def writeData(dat:pd.DataFrame) -> None:
    config = get_config()
    outfile = os.path.join(get_database_path(), config["databases"]["ModelSeed"]["reformat"])
    dat.to_csv(outfile)

def extractModelSEEDAnnotationsFromAlias(alias:str) -> tuple[dict,list]:
    '''
    Takes an entry from the model seed database ("https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/compounds.tsv") and retrieves the crosslinked databases into the annotation dictionary
    '''

    # the string is basically already a dictionary, it just needs to be passed back to it
    # first we need to replace all possible quotation marks
    if type(alias) == str:
        # there are some inconsistencies in the data which needs to be handled
        alias = alias.replace(" : "," ")
        alias = alias.replace("[|]","")
        alias = alias.replace('"',"''")
        new_alias = '''{"'''+alias.replace(''': ''','''":["''').replace('''|''','''"],"''').replace('''; ''','''","''')+'''"]}''' # the ''' are awkward, but necessary, because of ' sometimes denotes something like 3'-triphosphate, which makes ' as delimiter unusable
        new_anno = eval(new_alias)

        # get the alias names and remove them from the dictionary
        if "Name" in new_anno.keys():
            names = new_anno["Name"]
            del new_anno["Name"]
        else:
            names = list()

        # correct the keys to conform to the identifiers.org 
        new_anno = correctAnnotationKeys(new_anno)
    else:
        new_anno = dict()
        names = list()
    return(new_anno, names)

def correctAnnotationKeys(anno:dict, allow_missing_dbs: bool = False) -> dict:
    "Takes an annotation dictionary and tries to correct the keys in the dictionary to conform with the identifiers.org prefixes"

    # try to find the correct identifiers from identifiers.org
    ## first get possible prefixes
    prefixes = identifier_prefixes

    ## now check if the keys of the annos can be found in the prefixes
    new_anno_corrected = dict()
    for key, value in anno.items():
        # handle kegg drug vs. compound
        if key.lower() == "kegg":
            # sort the values into compounds and drugs
            drugs = []
            compounds = []
            for val in value:
                if val.startswith("D"):
                    drugs.append(val)
                elif val.startswith("C"):
                    compounds.append(val)
            if len(drugs) >0:
                new_anno_corrected["kegg.drug"] = drugs
            if len(compounds) >0:
                new_anno_corrected["kegg.compound"] = compounds
        # remove undefined pubchem entries - compounds and substances can have the same ID, but are usually different metabolites 
        if key.lower() == "pubchem":
            warnings.warn("Removing pubchem entries because of missing information about compound or substance:" + str(value))
        # handle database specific stuff like the .compound suffix
        elif key.lower()+".compound" in prefixes:
            new_anno_corrected[key.lower() + ".compound"] = value
        elif key.lower()+".metabolite" in prefixes:
            new_anno_corrected[key.lower() + ".metabolite"] = value
        elif key.lower()+".cpd" in prefixes:
            new_anno_corrected[key.lower() + ".cpd"] = value
        elif key.lower()+".chemical" in prefixes:
            new_anno_corrected[key.lower() + ".substance"] = value
        elif key.lower()+".substance" in prefixes:
            new_anno_corrected[key.lower() + ".chemical"] = value
        elif key.lower()+".entity" in prefixes:
            new_anno_corrected[key.lower() + ".entity"] = value
        elif key.lower()+".ligand" in prefixes:
            new_anno_corrected[key.lower() + ".ligand"] = value
        elif key.lower()+".smallmolecule" in prefixes:
            new_anno_corrected[key.lower() + ".ligand"] = value
        elif key.lower()+".inhibitor" in prefixes:
            new_anno_corrected[key.lower() + ".inhibitor"] = value
        elif key.lower()+".pthcmp" in prefixes:
            new_anno_corrected[key.lower() + ".pthcmp"] = value
        elif key.lower() in prefixes:
            new_anno_corrected[key.lower()] = value
    return(new_anno_corrected)

def handleNamesDBs(dat:pd.DataFrame) -> pd.DataFrame:
    # create a common name column from 'abbreviation','name','alias'
    n = dat.shape[0]
    name_db_dict = {"names" : [""]*n,
                    "DBs" : [""]*n}
    for i in range(n):
        # handle name column
        name = dat.loc[i,"name"]
        if type(name) != str:
            name = ""
        # handle abbreviation column
        abb = dat.loc[i,"abbreviation"]
        if type(abb) != str:
            abb = ""
        # handle aliases column
        databases,alt_names = extractModelSEEDAnnotationsFromAlias(dat.loc[i,"aliases"])

        # make all names unique, sort them and concatenate them as string seperated by "|"
        names = [name,abb]
        names.extend(alt_names)
        names = [x.lower() for x in names]
        names = list(set(names))
        names.sort()
        names = "|".join(names)

        # add names and database annotations
        name_db_dict["names"][i] = names
        name_db_dict["DBs"][i] = str(databases)

    # put everything in a dataframe
    name_db_df = pd.DataFrame(name_db_dict, index = dat.id)
    return(name_db_df)

def reformatModelSeed()->None:
    # a function to reformat the metabolite annotation table of ModelSEED to a standardized format
    # get the data
    dat = getData()
    # get names and databases
    namesDBs = handleNamesDBs(dat)
    # get inchis
    inchis = dat.smiles.apply(smile2inchi)
    inchis.name = "inchi"
    # combine everything back to a dataframe
    inchis.index = dat.id
    dat.index = dat.id
    dat_all = pd.concat([namesDBs,inchis,dat],axis =1)
    # write the data
    writeData(dat_all)
