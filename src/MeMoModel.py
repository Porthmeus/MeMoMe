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
from src.parseMetaboliteInfoFromSBML import parseMetaboliteInfoFromSBML, parseMetaboliteInfoFromSBMLMod



class MeMoModel():
    ''' The model class of MeMoMe. Core (for now) is a list of metabolites storing the relevant information. Further it will store the cobra representation for the model'''

    def __init__(self):
        ''' Stump for construction '''
        metabolites = []
        cobra_model = None
    
    @classmethod
    def fromPath(self, sbmlfile:Path) -> None:
        ''' Read the model from the SBML file '''
        self.metabolites = parseMetaboliteInfoFromSBML(sbmlfile, validate = True)


    @classmethod
    def fromModel(self, model:cb.Model) -> None:
        ''' Read the model from the cobra model '''
        raise ValueError("Parsing from cobrapy model not supported yet")

    @classmethod
    def fromSBML(self, model:sbml.model) -> None:
        ''' Read the model from the libsbml model '''
        self.metabolites = parseMetaboliteInfoFromSBMLMod(model)
    
