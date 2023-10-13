# Porthmeus
# 19.06.23

'''
Store relevant information of a SBML file in a MeMoModel
'''
from __future__ import annotations

from pathlib import Path

import cobra as cb
import libsbml as sbml
import pandas as pd
from src.annotateBulkRoutines import *
from src.matchInchi import matchInchi
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
        cobra_model = cb.io.read_sbml_model(str(sbmlfile))
        # TODO CATCH ERROR FROM COBRA
        _id = cobra_model.id
        return MeMoModel(metabolites=metabolites, cobra_model=cobra_model, _id=_id)

    @classmethod
    def fromModel(cls, model: cb.Model) -> MeMoModel:
        """ Read the model from the cobra model """
        cobra_model = model
        id = model.id
        metabolites = parseMetaboliteInfoFromCobra(model)
        return MeMoModel(metabolites=metabolites, _id=id)

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
        unannoted, annoted_by_chebi = annotateChEBI(self.metabolites)
        print(f'Out of {unannoted} metabolites that don\'t have an INCHI string, {annoted_by_chebi} were annotated by chebi')
        # GO BULK WISE ThORUGH BIGG AND VMH AND MODELSEED, try to extract as much as possible
        annotateLove(self.metabolites)


    def compare(self, model2: MeMoModel) -> pd.DataFrame:
        """ compares the metabolites of two models and returns a data frame with additional information """
        res = self.compareOnInchi(model2)
        # TODO add data base comparison
        return(res)


    def compareOnInchi(self, model2: MeMoModel) -> pd.DataFrame:
        # start with the comparison of inchi strings
        mod1_inchis = pd.DataFrame({"met_id" : [x.id for x in [y for y in self.metabolites]],
                "inchis" : [x._inchi_string for x in [y for y in self.metabolites]]})
        mod2_inchis = pd.DataFrame({"met_id" : [x.id for x in [y for y in model2.metabolites]],
                "inchis" : [x._inchi_string for x in [y for y in model2.metabolites]]})
        # go through the inchis of the second model and find the corresponding inchi in self
        matches = {"met_id1" : [],
                "met_id2" : [],
                "inchi_string":[],
                "charge_diff" : []}
        for i in range(len(mod1_inchis)):
            inchi1 = mod1_inchis.loc[i, "inchis"]
            if inchi1 != None:
                id1 = mod1_inchis.loc[i,"met_id"]
                for j in range(len(mod2_inchis)):
                    inchi2 = mod2_inchis.loc[j, "inchis"]
                    if inchi2 != None:
                        id2 = mod2_inchis.loc[j,"met_id"]
                        res = matchInchi(inchi1, inchi2)
                        if res[0] == True:
                            matches["met_id1"].append(id1)
                            matches["met_id2"].append(id2)
                            matches["inchi_string"].append(inchi1)
                            matches["charge_diff"].append(res[1])
        inchiRes = pd.DataFrame(matches)
        return(inchiRes)


