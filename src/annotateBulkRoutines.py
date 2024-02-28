# Porthmeus
# 07.07.23

'''
This file contains function to annotated MeMoMetabolites in a bulk manner. This means a list of MeMoMetabolites is parsed and the annotation takes place on the metabolites in that list. This is mainly to use within the MeMoModel.annotate() function
'''

import os
import numpy as np
import pandas as pd
import json
import warnings
from pathlib import Path
from src.download_db import get_config, get_database_path
from src.annotateInchiRoutines import smile2inchi
from io import StringIO
import re

from src.MeMoMetabolite import MeMoMetabolite
from src.annotateInchiRoutines import findOptimalInchi
from src.download_db import databases_available, get_config, get_database_path
from src.parseMetaboliteInfos import getAnnoFromIdentifierURL

def annotateChEBI(metabolites: list[MeMoMetabolite]) -> tuple[int, int]:
    """ Annotate the metaboltes with Inchis from ChEBI """

    # check if databases are available and if so load the config
    if databases_available():
        config = get_config()
        db_path = get_database_path()

    # check if any unannotated metabolites exist and whether these have a ChEBI entry
    ids = [x for x, y in enumerate(metabolites) if y._inchi_string == None and "chebi" in y.annotations.keys()]
    not_annotated_metabolites = len(ids)
    annotated_by_chebi = 0
    # check if chebis ids are actually present in the annotation slot
    # TODO This does not work with 35, because we try to access the "chebi" even for thos metabolies that don;t have one
    annos = any(["chebi" in x.annotations.keys() for i, x in enumerate(metabolites) if i in ids])
    if annos:
        # download the information from the server
        #chebi_db = pd.read_table("/home/td/Downloads/chebiId_inchi.tsv")
        chebi_db = pd.read_table(db_path.joinpath(config["databases"]["ChEBI"]["file"]))
        chebi_db.index = chebi_db['CHEBI_ID']

        # annotate the metabolites with the inchi_string
        for i in ids:
            # For each metabolite that does not have a INCHI string get its chebi id #TODO  KEY OR STRING
            chebis = [int(x.replace("CHEBI:", "")) for x in metabolites[i].annotations["chebi"]]
            # Find all the corresponding INCHI key
            chebis = [x for x in chebis if x in chebi_db.index]
            inchis = np.unique(chebi_db.loc[chebis, "InChI"])
            if len(inchis) > 0:
                inchi = findOptimalInchi(inchis)
                annotated_by_chebi += 1

            metabolites[i].set_inchi_string(inchi)

    return not_annotated_metabolites, annotated_by_chebi


def annotateLove(metabolites: list[MeMoMetabolite]) -> tuple[int, int]:
    """ Annotate the metaboltes with Inchis from ChEBI """

    # check if any unannotated metabolites exist
    ids = [x for x, y in enumerate(metabolites) if y._inchi_string == None]
    not_annotated_metabolites = len(ids)
    print("NOT ANNO", not_annotated_metabolites)

    return 0, 0

def annotateBiGG_entry(entry:str, annotations:dict = dict(), database:pd.DataFrame = pd.DataFrame()) -> dict:
    """
    A small helper function to avoid redundant code
    Uses a BiGG identifiers and annotates it with the identifiers.org entries.
    The results are stored int he annotations dictionary. This is empty by
    default, but if it exists, the entries will be simply added (duplicates are
    allowed here).
    database - is either empty (default) or a pandas data.frame with the data from the BiGG model homepage for the metabolites. If it is empty, the function will try to load the database from the config file
    """
    # check if the database was given, if not, try to load it
    if len(database) == 0:
        # load the database
        config = get_config()
        db_path =  os.path.join(get_database_path(), config["databases"]["BiGG"]["file"])
        bigg = pd.read_table(db_path)
    else:
        bigg = database

    # extract the relevant annotation information and return it
    urls = bigg.loc[bigg["universal_bigg_id"]==entry,"database_links"]
    # check if there are entries which are not NA
    urls = urls.loc[~pd.isna(urls)]
    if len(urls) >0:
        urls =";".join(list(urls))
        urls = urls.split(";")
        for url in urls:
            src,val =getAnnoFromIdentifierURL(url)
            if src in annotations.keys():
                annotations[src].append(val)
            else:
                annotations[src] = [val]
    
    # extract the saved name in the database and return it
    names = list(pd.unique(bigg.loc[bigg["universal_bigg_id"]==entry,"name"]))
    return(annotations, names)

    
    

def annotateBiGG(metabolites: list[MeMoMetabolite]) -> tuple[int, int]:
    """
    Annotate a list of metabolites with the entries from BiGG. Look for BiGG ids in the annotation slot of the metabolites and if one is found use these. Since BiGG does not provide any InChI strings, this will not results in any new annotated Inchi strings, but it will increase the number of entries in the metabolite annotation.
    """
    
    # load the database
    config = get_config()
    db_path =  os.path.join(get_database_path(), config["databases"]["BiGG"]["file"])
    bigg = pd.read_table(db_path)
    
    counter = 0
    # go through the metabolites and check if there is data which can be added
    for met in metabolites:
        new_met_anno = dict()
        if "bigg.metabolite" in met.annotations.keys():
            for entry in met.annotations["bigg.metabolite"]:
                new_met_anno, new_names = annotateBiGG_entry(entry = entry,
                        annotations = new_met_anno,
                        database = bigg)
                # add new names
                met.add_names(new_names)

            # remove duplicates from the annotations
            for key in new_met_anno.keys():
               new_val = list(set(new_met_anno[key]))
               new_met_anno[key] = new_val
            
            # add the annotations to the slot in the metabolites
            if len(new_met_anno) > 0:
                met.add_annotations(new_met_anno)
                counter = counter +1
    
    # return a message how many annotations could be found
    print("Found " + str(counter) + " annotations in the BiGG database, which have been added")

    # get the tuple for returning the annotation counts - this is just done to comply to a standard (there will be no InChI annotation from this database
    ret1,ret2 = annotateLove(metabolites)
    return ret1,ret2


def annotateBiGG_id(metabolites: list[MeMoMetabolite]) -> tuple[int, int]:
    """
    Annotate a list of metabolites with the entries from BiGG. Look for BiGG ids in the metabolite._id slot and if one is found use these. Since BiGG does not provide any InChI strings, this will not results in any new annotated Inchi strings, but it will increase the number of entries in the metabolite annotation.
    """

    # load the database
    config = get_config()
    db_path =  os.path.join(get_database_path(), config["databases"]["BiGG"]["file"])
    bigg = pd.read_table(db_path)
    
    counter = 0
    for met in metabolites:
        new_met_anno = dict()
        #if met._id in bigg["universal_bigg_id"]: # does not work, for whatever reason
        if any(bigg["universal_bigg_id"]==met._id):
            new_met_anno,new_names = annotateBiGG_entry(entry = met._id,
                    annotations = new_met_anno,
                    database = bigg)
            # add new names
            met.add_names(new_names)

            # remove duplicates from the annotations
            for key in new_met_anno.keys():
               new_val = list(set(new_met_anno[key]))
               new_met_anno[key] = new_val
        
            # add the annotations to the slot in the metabolites
            if len(new_met_anno) > 0:
                met.add_annotations(new_met_anno)
                counter = counter +1

    # return a message how many annotations could be found
    print("Found " + str(counter) + " annotations in the BiGG database, which have been added")

    # get the tuple for returning the annotation counts - this is just done to comply to a standard (there will be no InChI annotation from this database
    ret1,ret2 = annotateLove(metabolites)
    return ret1,ret2

def extractModelSeedAnnotationsFromAlias(alias:str) -> tuple[dict,list]:
    '''
    Takes an entry from the model seed database ("https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/compounds.tsv") and retrieves the crosslinked databases into the annotation dictionary
    '''

    # the string is basically already a dictionary, it just needs to be passed back to it
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
    return(new_anno, names)

def correctAnnotationKeys(anno:dict) -> dict:
    "Takes an annotation dictionary and tries to correct the keys in the dictionary to conform with the identifiers.org prefixes"

    ## TODO: currently the data base is read everytime we run this function (which could be several k times) - find a better solution for that
    # try to find the correct identifiers from identifiers.org
    ## first get possible prefixes

    config = get_config()
    identifiers_file = os.path.join(get_database_path(), config["databases"]["Identifiers"]["file"])
    with open(identifiers_file, "r") as json_file:
        identifiers = json.load(json_file)
    prefixes = json.dumps(identifiers["_embedded"]["namespaces"])
    prefixes = pd.read_json(StringIO(prefixes))
    prefixes = list(prefixes["prefix"])

    ## now check if the keys of the annos can be found in the prefixes
    new_anno_corrected = dict()
    for key, value in anno.items():
        if key.lower() in prefixes:
            new_anno_corrected[key.lower()] = value
        # handle kegg drug vs. compound
        elif key.lower() == "kegg":
            # sort the values into compounds and drugs
            drugs = []
            compounds = []
            for val in value:
                if val.startswith("D"):
                    drugs.append[val]
                elif val.startswith("C"):
                    compounds.append[val]
            if len(drugs) >0:
                new_anno_corrected["kegg.drug"] = drugs
            if len(compounds) >0:
                new_anno_corrected["kegg.compound"] = compounds
        # remove undefined pubchem entries - compounds and substances can have the same ID, but are usually different metabolites 
        elif key.lower() == "pubchem":
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
        
    return(new_anno_corrected)

def extractModelSeedpKapKb(pkab:str) -> dict:
    """Auxiliary function to reformat the values in the modelseed database for pKa and pKb into dictionaries"""
    # make sure there is an actual value for the entry
    if not pd.isna(pkab):
        # split the string into individual pkA/Bs
        pkab_split = pkab.split(";")

        # sort the coordinates and the pkA/B value into a dictonary
        pkab_dict = {}
        for entry in pkab_split:
            coord = re.sub(r"^(.*):(.*):(.*)$", r"\g<1>:\g<2>",entry)
            pkab_val = float(re.sub(r"^(.*):(.*):(.*)$", r"\g<3>",entry))
            pkab_dict[coord] = pkab_val
    else:
        pkab_dict={"0:0",float("nan")}

    return(pkab_dict)


def annotateModelSeed_entry(entry:str,  database:pd.DataFrame = pd.DataFrame()) -> list[dict, list, dict, dict]:
    """
    A small helper function to avoid redundant code
    Uses a ModelSeed identifiers and annotates it with the identifiers.org entries.
    The results are stored in the annotations dictionary. This is empty by
    default, but if it exists, the entries will be simply added (duplicates are
    allowed here). Further it will return a list of alternative names, and pKa/pKb values
    database - is either empty (default) or a pandas data.frame with the data from the modelseed homepage for the metabolites. If it is empty, the function will try to load the database from the config file
    """
    
    # check if the database was given, if not, try to load it
    if len(database) == 0:
        # load the database
        config = get_config()
        db_path =  os.path.join(get_database_path(), config["databases"]["ModelSeed"]["file"])
        mseed = pd.read_table(db_path)
    else:
        mseed = database

    # extract the relevant annotation information and return it
    aliases = mseed.loc[mseed["id"]==entry,"aliases"]
    # check if there are entries which are not NA
    aliases = aliases.loc[~pd.isna(aliases)]
    # get the information from the database
    new_anno = dict()
    new_names = list()
    if len(aliases) >0:
        for alias in aliases:
            annos, names = extractModelSeedAnnotationsFromAlias(alias)
            # combine annos for the entry
            for key,value in annos.items():
                if key in new_anno.keys():
                    new_anno[key].extend(value)
                else:
                    new_anno[key] = value
            # combine names
            new_names.extend(names)
        
        # remove duplicates
        for key in new_anno.keys():
            new_anno[key] = list(set(new_anno[key]))
        new_names = list(set(new_names))
    
    # additionally extract the pKa and pKb values
    pka = mseed.loc[mseed["id"]==entry,"pka"]
    pkb = mseed.loc[mseed["id"]==entry,"pkb"]
    # check if there are entries which are not NA
    pka = pka.loc[~pd.isna(pka)]
    pkb = pkb.loc[~pd.isna(pkb)]
    # get the pka values
    if len(pka) > 0:
        pka_vals = dict()
        for val in pka:
            pka_temp = extractModelSeedpKapKb(val)
            for key,value in pka_temp.items():
                if key in pka_vals.keys():
                    pka_vals[key] = (pka_vals[key] + value)/2
                else:
                    pka_vals[key] = value
    else:
        pka_vals = {}
    # do the same for pkb
    if len(pkb) > 0:
        pkb_vals = dict()
        for val in pkb:
            pkb_temp = extractModelSeedpKapKb(val)
            for key,value in pkb_temp.items():
                if key in pkb_vals.keys():
                    pkb_vals[key] = (pkb_vals[key] + value)/2
                else:
                    pkb_vals[key] = value
    else:
        pkb_vals = {}

    return(new_anno, new_names, pka_vals, pkb_vals)

def annotateModelSeed_id(metabolites: list[MeMoMetabolite]) ->int:
    """
    Annotate a list of metabolites with the entries from ModelSeed. Look for ModelSeed ids in the metabolite._id slot and if one is found use these. Since ModelSeed does provide any InChI strings, pKs and pkB, all these elements will be added to the MeMoMetabolites, if possible"""


    # load the database
    config = get_config()
    db_path =  os.path.join(get_database_path(), config["databases"]["ModelSeed"]["file"])
    mseed = pd.read_table(db_path)
    
    counter = [0,0,0,0,0] # counter for names, annotation, inchi_string, pka, pkb
    for met in metabolites:
        if any(mseed["id"]==met._id):
            
            # get the names, annotations, pka and pkb
            new_met_anno,new_names,new_pka,new_pkb = annotateModelSeed_entry(entry = met._id,
                    annotations = new_met_anno,
                    database = bigg)

            # get the inchi strings
            smiles = mseeds.loc[mseed["id"] == met._id, smiles]
            smiles = smiles.loc[~pd.isna(smiles)]
            inchi_strings = []
            # there are only smiles in modelseed and rdkit sometimes fails to convert them - report here if that happens
            for smile in smiles:
                try:
                    inchi_strings.append(smile2inchi(smile))
                except Exception as e:
                    warnings.warn("Could not convert smile to inchi for metabolite {met} and {smile}\n". format(met = met._id, smile = smile) + str(e))
            # get the correct inchi_string, if there was more than one
            if len(inchi_strings) > 1:
                inchi_string = findOptimalInchi(inchi_strings)
            elif len(inchi_strings) == 1:
                inchi_string = inchi_strings[0]
            else:
                inchi_string = None

            # add new names
            if len(new_names) > 0:
                #counter[0] = counter[0] + 
                met.add_names(new_names)

            # add the annotations to the slot in the metabolites
            if len(new_met_anno) > 0:
                #counter[1] = counter[1] + 
                met.add_annotations(new_met_anno)
            
            # add the inchi_string
            if inchi_string != None:
                counter[2] = counter[2] + met.add_inchi_string(inchi_string)

            # add the pka
            if len(new_pka) > 0:
                # counter[3] = counter[3] +
                met.add_pka(new_pka)

            # add the pkb
            if len(new_pkb) > 0:
                # counter[4] = counter[4] + 
                met.add_pkb(new_pkb)

    return counter[2]

def annotateModelSeed(metabolites: list[MeMoMetabolite]) ->int:
    """
    Annotate a list of metabolites with the entries from ModelSeed. Look for ModelSeed ids in the metabolite.anotation slot and if one is found use these. Since ModelSeed does provide any InChI strings, pKs and pkB, all these elements will be added to the MeMoMetabolites, if possible"""


    # load the database
    config = get_config()
    db_path =  os.path.join(get_database_path(), config["databases"]["ModelSeed"]["file"])
    mseed = pd.read_table(db_path)
    
    counter = [0,0,0,0,0] # counter for names, annotation, inchi_string, pka, pkb
    for met in metabolites:
        if "seed.compound" in met.annotation.keys():
            mseeds = met.annotations["seed.compounds"]
            met_counter = [0,0,0,0,0]
            for seed_id in mseeds:
                if any(mseed["id"]==seed_id):
                    # get the names, annotations, pka and pkb
                    new_met_anno,new_names,new_pka,new_pkb = annotateModelSeed_entry(
                            entry = met._id,
                            annotations = new_met_anno,
                            database = bigg)

                    # get the inchi strings
                    smiles = mseeds.loc[mseed["id"] == met._id, smiles]
                    smiles = smiles.loc[~pd.isna(smiles)]
                    inchi_strings = []
                    # there are only smiles in modelseed and rdkit sometimes fails to convert them - report here if that happens
                    for smile in smiles:
                        try:
                            inchi_strings.append(smile2inchi(smile))
                        except Exception as e:
                            warnings.warn("Could not convert smile to inchi for metabolite {met} and {smile}\n". format(met = met._id, smile = smile) + str(e))
                    # get the correct inchi_string, if there was more than one
                    if len(inchi_strings) > 1:
                        inchi_string = findOptimalInchi(inchi_strings)
                    elif len(inchi_strings) == 1:
                        inchi_string = inchi_strings[0]
                    else:
                        inchi_string = None
                

                    # add new names
                    if len(new_names) > 0:
                        #met_counter[0] = met_counter[0] + 
                        met.add_names(new_names)

                    # add the annotations to the slot in the metabolites
                    if len(new_met_anno) > 0:
                        #met_counter[1] = met_counter[1] + 
                        met.add_annotations(new_met_anno)
                    
                    # add the inchi_string
                    if inchi_string != None:
                        met_counter[2] = met_counter[2] + met.add_inchi_string(inchi_string)

                    # add the pka
                    if len(new_pka) > 0:
                        # met_counter[3] = met_counter[3] +
                        met.add_pka(new_pka)

                    # add the pkb
                    if len(new_pkb) > 0:
                        # met_counter[4] = met_counter[4] + 
                        met.add_pkb(new_pkb)

            # sum the counters
            met_counter = [int(x>0) for x in met_counter]
            counter = [x+y for x,y in zip(counter, met_counter)]

    return counter[2]
