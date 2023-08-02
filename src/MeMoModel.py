# Porthmeus
# 19.06.23

'''
Store relevant information of a SBML file in a MeMoModel
'''
from __future__ import annotations

from pathlib import Path

import cobra as cb
import libsbml as sbml
from src.annotateBulkRoutines import *
from src.parseMetaboliteInfos import parseMetaboliteInfoFromSBML, parseMetaboliteInfoFromSBMLMod, \
    parseMetaboliteInfoFromCobra


class MeMoModel:
    """ The model class of MeMoMe. Core (for now) is a list of metabolites storing the relevant information. Further
    it will store the cobra representation for the model"""

    def __init__(self,
                 metabolites: list[MeMoMetabolite] = [],
                 cobra_model: cb.Model = None,
                 _id: str = None) -> None:
        ''' Stump for construction '''
        self.metabolites = metabolites
        self.cobra_model = cobra_model
        self._id = _id

    @classmethod
    def fromPath(cls, sbmlfile: Path) -> MeMoModel:
        """ Read the model from the SBML file """
        metabolites = parseMetaboliteInfoFromSBML(sbmlfile, validate=True)
        cobra_model = cb.io.read_sbml_model(sbmlfile)
        _id = cobra_model.id
        return MeMoModel(metabolites=metabolites, cobra_model=cobra_model, _id=_id)

    @classmethod
    def fromModel(cls, model: cb.Model) -> MeMoModel:
        """ Read the model from the cobra model """
        raise NotImplementedError("Parsing from cobrapy model not supported yet")
        # cobra_model = model
        # id = model.id
        # metabolites = parseMetaboliteInfoFromCobra(model)
        # return MeMoModel(metabolites=metabolites, _id=id)

    @classmethod
    def fromSBML(cls, model: sbml.Model) -> MeMoModel:
        """ Read the model from the libsbml model """
        metabolites = parseMetaboliteInfoFromSBMLMod(model)
        _id = model.getId()
        return MeMoModel(metabolites=metabolites, _id=_id)

    def annotate(self) -> None:
        """Goes through the different bulk annotation methods and tries to annotate InChI strings to the metabolites
        in the model"""
        # Use ChEBI
        annotateChEBI(self.metabolites)
