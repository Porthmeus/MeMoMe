# Porthmeus
# 15.02.23

'''
Read an SBML file and return a pandas df with all the information on its metabolites
'''

import libsbml as sbml
import re
import pandas as pd
import warnings
from os.path import exists
from pathlib import Path
from src.MeMoMetabolite import MeMoMetabolite

def getAnnotationFromMet(met:sbml.Species) -> dict:
    ''' extract a dictionary containing additional annotations'''
    results = dict()
    urls = met.getAnnotationString().split("\n")
    if len(urls) > 1:
        anno_source = [x.split("/")[3].strip('"') for x in urls if re.search("identifiers.org",x) != None]
        anno_val = ["/".join(x.split("/")[4:]).replace('"/>',"") for x in urls if re.search("identifiers.org",x) != None]
        for src_i in range(len(anno_source)):
            src = anno_source[src_i]
            val = anno_val[src_i]
            if src in results.keys():
                results[src].append(val)
            else:
                results[src] = [val]
    return(results)
                
def validateSBML(modelDoc:sbml.SBMLDocument, strict:bool = False) -> None:
    '''
    Uses the libsbml validator to validate an sbml file. Throws an error or warning if the validation fails.
    Params:
        modelDoc:libsbml.SBMLDoc - a Path object to the sbml file.
        strict:bool - Whether to raise an error or just a warning. Default: False - only raise Warning.
    '''
    validator = sbml.SBMLValidator()
    errors = validator.validate(modelDoc)
    if errors > 0:
        text = '''{MSG}: SBML file not valid, following failures have been found:\n''' + validator.getErrorLog().printErrors()
        if strict:
            raise IOError(text.format(MSG = "ERROR"))
        else:
            warning.warn(text.format(MSG = "WARNING"))
    
def annotateMetabolites_bulk(metabolites:list[MeMoMetabolite]) -> None:
    '''
    This function tries to get the InChI string for each metabolite by leveraging the information in the annotation attribute of the metabolite
    '''


def parseMetaboliteInfoFromSBML(modelfile: Path, validate:bool = True) -> list:
    '''
    Parses through the SBML and extracts metabolite objects from it. 
    Params:
        modelfileTakes:Path - a Path object to a valid sbml file
        validate:bool - whether to validate the modelfile before parsing, default: True.

    '''
    # handle non-existing files as sbml won't complain
    if not modelfile.exists():
        raise FileExistsError(str(modelfile) + " does not exist")

    # read the sbml file
    reader = sbml.SBMLReader()
    doc = reader.readSBMLFromFile(modelfile)

    # validate if necessary
    if validate:
        validateSBML(doc, strict = False)

    # get the actual model and its id from the file
    mod = doc.getModel()
    mod_id = mod.getId()

    # go through the metabolites in the model and extract the relevant information
    metabolites = mod.getListOfSpecies() 
    memoMets = []
    for met in metabolites:
        annotations = getAnnotationFromMet(met)
        # check if the inchi_string is available
        if "inchi" in annotations.keys():
            inchi_string = annotations['inchi'][0]
        else:
            inchi_string = None
        # create a MeMoMetabolite
        memomet = MeMoMetabolite(_id = met.getId(),
                orig_ids = [met.getId()],
                _model_id = mod_id,
                names = [met.getName()],
                _inchi_string = inchi_string,
                _formula = None,
                _charge = met.getCharge(),
                annotations = annotations)

        # check if the MeMoMetabolite is already in the list of metabolites and whether it can be merged.
        if memomet._id in [x._id for x in memoMets]:
            idx = [i for i in range(len(memoMets)) if memoMets[i]._id == memomet._id][0]
            memoMets[idx].merge(memomet)
        else:
            memoMets.append(memomet)

    return(memoMets)