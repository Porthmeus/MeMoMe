# Porthmeus
# 15.02.23
# edit: 07.07.23

'''
Functions to parse relevant metabolite infos into a list of MeMoMetabolites from different objects like: a path to an SBML file, a libsbml.Model object or a cobra.Model object
'''

import libsbml as sbml
import re
import pandas as pd
import warnings
import numpy as np
import cobra as cb
from os.path import exists
from pathlib import Path
from rdkit import Chem

from debugpy.common.log import warning

from src.MeMoMetabolite import MeMoMetabolite
from src.annotateInchiRoutines import validateInchi, findOptimalInchi

def getAnnoFromIdentifierURL(url:str) -> tuple[str, str]:
    """Small function to extract the db and db id from an identifier.org url"""
    db = url.split("/")[3].strip('"')
    db_id = "/".join(url.split("/")[4:]).replace('"/>', "")
    return db, db_id

def getAnnotationFromMet(met: sbml.Species) -> dict:
    ''' extract a dictionary containing additional annotations'''
    results = dict()
    urls = met.getAnnotationString().split("\n")
    # filter those anntations which have an identifier.org url attached
    urls = [x for x in urls if re.search("identifiers.org", x) != None]
    # if annotations are available extract the identifier.org db and id and put it in the annotations dictionary
    if len(urls) > 0:
        for url in urls:
            src,val =getAnnoFromIdentifierURL(url)
            if src in results.keys():
                results[src].append(val)
            else:
                results[src] = [val]

    # remove duplicates from the annotations
    for key in results.keys():
       new_val = list(set(results[key]))
       results[key] = new_val
    return results


def validateSBML(modelDoc: sbml.SBMLDocument, strict: bool = False) -> None:
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
            raise IOError(text.format(MSG="ERROR"))
        else:
            warning.warn(text.format(MSG="WARNING"))


def parseMetaboliteInfoFromSBML(modelfile: Path, validate: bool = True) -> list:
    '''
    Parses through the SBML and extracts metabolite objects from it. 
    Parameters
    ----------
        modelfile : Path - a Path object to a valid sbml file
        validate : bool - whether to validate the modelfile before parsing, default: True.

    '''
    # handle non-existing files as sbml won't complain
    if not modelfile.exists():
        raise FileExistsError(str(modelfile) + " does not exist")

    # read the sbml file
    reader = sbml.SBMLReader()
    doc = reader.readSBMLFromFile(modelfile)

    # validate if necessary
    if validate:
        validateSBML(doc, strict=False)

    # get the actual model and its id from the file
    mod = doc.getModel()
    mod_id = mod.getId()

    # extract the metabolites 
    memoMets = parseMetaboliteInfoFromSBMLMod(model=mod)
    return memoMets


def parseMetaboliteInfoFromSBMLMod(model: sbml.Model) -> list:
    '''
    Parses through the SBML model object and extracts metabolite objects from it. 
    Parameters
    ----------
        model : sbml.model -  a valid sbml model which can be parsed
    '''

    # go through the metabolites in the model and extract the relevant information
    metabolites = model.getListOfSpecies()
    memoMets = []
    for met in metabolites:
        annotations = getAnnotationFromMet(met)
        # check if the inchi_string is available
        if "inchi" in annotations.keys():
            inchi_string = annotations['inchi'][0]
            # check for validity of the inchi
            if not validateInchi(inchi_string):
                inchi_string = None
        else:
            inchi_string = None
        # create a MeMoMetabolite
        memomet = MeMoMetabolite(_id=met.getId(),
                                 orig_ids=[met.getId()],
                                 _model_id=model.getId(),
                                 names=[met.getName()],
                                 _inchi_string=inchi_string,
                                 _formula=None,
                                 _charge=met.getCharge(),
                                 annotations=annotations)

        # check if the MeMoMetabolite is already in the list of metabolites and whether it can be merged.
        if memomet.id in [x.id for x in memoMets]:
            idx = [i for i in range(len(memoMets)) if memoMets[i]._id == memomet._id][0]
            memoMets[idx].merge(memomet)
        else:
            memoMets.append(memomet)

    return memoMets


def parseMetaboliteInfoFromCobra(model: cb.Model) -> list[MeMoMetabolite]:
    '''
    Parses through and cobra model object and extracts metabolite objects from it. 
    Parameters
    ----------
        model : cobra.model -  a valid sbml model which can be parsed
    '''

    metabolites = []
    memoMets = []
    for met in model.metabolites:
        annotations = met.annotation
        # make lists out of the annotations
        for key in annotations.keys():
            if type(annotations[key]) != list:
                annotations[key] = [annotations[key]]

        if "inchi" in annotations.keys():
            inchi = annotations["inchi"]
            # check whether there is more than one inchi in the annotations

            if len(inchi) > 1:
                inchi = findOptimalInchi(inchi, charge = met.charge)
                if inchi is None:
                    raise NotImplementedError()
            elif validateInchi(inchi[0]):
                inchi = inchi[0]
            else:
                inchi = None
        else:
            inchi = None
        # create the MeMoMetabolite object
        memomet = MeMoMetabolite(_id=met.id,
                                 orig_ids=[met.id],
                                 _model_id=[model.id],
                                 names=[met.name],
                                 _inchi_string=inchi,
                                 _formula=met.formula,
                                 _charge=met.charge,
                                 annotations=annotations)

        # check if the MeMoMetabolite is already in the list of metabolites and whether it can be merged.
        if memomet._id in [x._id for x in memoMets]:
            idx = [i for i in range(len(memoMets)) if memoMets[i]._id == memomet._id][0]
            memoMets[idx].merge(memomet)
        else:
            memoMets.append(memomet)

    return (memoMets)
