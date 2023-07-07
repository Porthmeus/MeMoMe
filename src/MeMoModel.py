# Porthmeus
# 19.06.23

'''
Store relevant information of a SBML file in a MeMoModel
'''



import libsbml as sbml
import re
import pandas as pd
import warnings
import numpy as np
import cobra as cb
from os.path import exists
from pathlib import Path
from src.MeMoMetabolite import MeMoMetabolite
from src.annotateInchiRoutines import findOptimalInchi
from src.parseMetaboliteInfos import parseMetaboliteInfoFromSBML, parseMetaboliteInfoFromSBMLMod, parseMetaboliteInfoFromCobra
from src.annotateBulkRoutines import *



class MeMoModel():
    ''' The model class of MeMoMe. Core (for now) is a list of metabolites storing the relevant information. Further it will store the cobra representation for the model'''

    def __init__(self,
            metabolites:list[MeMoMetabolite] = [],
            cobra_model:cb.Model = None,
            _id:str = None) -> None:
        ''' Stump for construction '''
        self.metabolites = metabolites
        self.cobra_model = cobra_model
        self._id = _id
    
    @classmethod
    def fromPath(self, sbmlfile:Path) -> None:
        ''' Read the model from the SBML file '''
        self.metabolites = parseMetaboliteInfoFromSBML(sbmlfile, validate = True)
        self.cobra_model = cb.io.read_sbml_model(sbmlfile)
        self._id = self.cobra_model.id

    @classmethod
    def fromModel(self, model:cb.Model) -> None:
        ''' Read the model from the cobra model '''
        self.cobra_model = model
        self._id = model.id
        self.metabolites = parseMetaboliteInfoFromCobra(model)
        raise ValueError("Parsing from cobrapy model not supported yet")

    @classmethod
    def fromSBML(self, model:sbml.Model) -> None:
        ''' Read the model from the libsbml model '''
        self.metabolites = parseMetaboliteInfoFromSBMLMod(model)
        self._id = model.getId()
    
    def annotate(self) -> None:
        '''Goes through the different bulk annotation methods and tries to annotate InChI strings to the metabolites in the model'''
        # Use ChEBI
        annotateChEBI(self.metabolites)

