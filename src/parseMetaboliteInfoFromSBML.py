# Porthmeus
# 15.02.23

'''
Read an SBML file and return a pandas df with all the information on its metabolites
'''


import libsbml as sbml
import re
import pandas as pd
from os.path import exists
from pathlib import Path

def parseMetaboliteInfoFromSBML(modelfile: Path) -> pd.DataFrame:
    '''
    This function reads an SBML file an returns a data frame with relevant information on all metabolites in the model
    '''
    # handle non-existing files as sbml won't complain
    if not modelfile.exists():
        raise FileExistsError(str(modelfile) + " does not exist")
    # read the sbml file
    reader = sbml.SBMLReader()
    doc = reader.readSBMLFromFile(modelfile)
    #doc = reader.readSBMLFromFile("./Recon3DModel_301.xml")
    mod = doc.getModel()
    # go through the metabolites in the model and extract the relevant information
    metabolites = mod.getListOfSpecies() 

    met_dic = {"ID":[],
            "Name":[],
            "Charge":[],
            "Compartment":[],
            "Chemical_formula":[],
            "Annotation_source":[],
            "Annotation_value":[]}
    for met_i in range(len(metabolites)):
        # get annotations first, if there are more than 1, use the length as multiplier
        met = mod.getSpecies(met_i)
        urls = met.getAnnotationString().split("\n")
        if len(urls) > 0:
            anno_source = [x.split("/")[3].strip('"') for x in urls if re.search("identifiers.org",x) != None]
            met_dic["Annotation_source"].extend(anno_source)
            met_dic["Annotation_value"].extend( ["/".join(x.split("/")[4:]).replace('"/>',"") for x in urls if re.search("identifiers.org",x) != None])
            mult = len(anno_source)
        else:
            mult = 1
            met_dic["Annotation_source"].extend([""] * mult)
            met_dic["Annotation_value"].extend( [""]* mult)
        
        # get the rest of the information
        met_dic["ID"].extend([met.getId()] * mult)
        met_dic["Name"].extend([met.getName()] * mult)
        met_dic["Charge"].extend([met.getCharge()] * mult)
        met_dic["Compartment"].extend([met.getCompartment()] * mult)
        
        # finally look for the chemical formula
        fbc = met.getPlugin("fbc")
        met_dic["Chemical_formula"].extend([fbc.getChemicalFormula()] * mult)

        # debugging
       # nrows = [len(x) for x in met_dic.values()]
       # if not all(x == nrows[0] for x in nrows):
       #     raise 
       #     break
    
    # cast the dictionary into a data frame and return it
    df = pd.DataFrame(met_dic)
    return(df)

def parseReactionInfoFromSBML(modelfile: Path) -> pd.DataFrame:
    '''
    This function reads an SBML file an returns a data frame with relevant information on all reaction in the model
    '''
    # handle non-existing files as sbml won't complain
    if not modelfile.exists():
        raise FileExistsError(str(modelfile) + " does not exist")
    # read the sbml file
    reader = sbml.SBMLReader()
    doc = reader.readSBMLFromFile(modelfile)
    #doc = reader.readSBMLFromFile("./Recon3DModel_301.xml")
    mod = doc.getModel()
    # go through the reactions in the model and extract the relevant information
    reactions = mod.getListOfReactions() 

    rxn_dic = {"ID":[],
            "Name":[],
            "Annotation_source":[],
            "Annotation_value":[]}
    for rxn_i in range(len(reactions)):
        # get annotations first, if there are more than 1, use the length as multiplier
        rxn = mod.getReaction(rxn_i)
        urls = rxn.getAnnotationString().split("\n")
        if len(urls) > 0:
            anno_source = [x.split("/")[3].strip('"') for x in urls if re.search("identifiers.org",x) != None]
            rxn_dic["Annotation_source"].extend(anno_source)
            rxn_dic["Annotation_value"].extend( ["/".join(x.split("/")[4:]).replace('"/>',"") for x in urls if re.search("identifiers.org",x) != None])
            mult = len(anno_source)
        else:
            mult = 1
            rxn_dic["Annotation_source"].extend([""] * mult)
            rxn_dic["Annotation_value"].extend( [""]* mult)
        
        # get the rest of the information
        rxn_dic["ID"].extend([rxn.getId()] * mult)
        rxn_dic["Name"].extend([rxn.getName()] * mult)
        

        # debugging
       # nrows = [len(x) for x in rxn_dic.values()]
       # if not all(x == nrows[0] for x in nrows):
       #     raise 
       #     break
    
    # cast the dictionary into a data frame and return it
    df = pd.DataFrame(rxn_dic)
    return(df)
