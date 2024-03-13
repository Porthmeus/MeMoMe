# Porthmeus
# 19.06.23

'''
Store relevant information of a SBML file in a MeMoModel
'''
from __future__ import annotations

from pathlib import Path
from warnings import warn
import cobra as cb
import libsbml as sbml
import pandas as pd
from src.annotateChEBI import annotateChEBI
from src.annotateBiGG import annotateBiGG, annotateBiGG_id
from src.annotateModelSEED import annotateModelSEED, annotateModelSEED_id
from src.annotateAux import AnnotationResult
from src.matchMets import matchMetsByDB, matchMetsByInchi, matchMetsByName
from src.parseMetaboliteInfos import parseMetaboliteInfoFromSBML, parseMetaboliteInfoFromSBMLMod, \
    parseMetaboliteInfoFromCobra
from src.MeMoMetabolite import MeMoMetabolite
from src.annotateInchiRoutines import inchiToMol, molToRDK, molToNormalizedInchi
from src.origin_databases import origin_databases
from rdkit import Chem
import logging

logger = logging.getLogger('logger')


class MeMoModel:
    """ The model class of MeMoMe. Core (for now) is a list of metabolites storing the relevant information. Further
    it will store the cobra representation for the model"""

    def __init__(self,
                 metabolites: list[MeMoMetabolite] = [],
                 cobra_model: cb.Model | None = None,
                 _id: str | None = None) -> None:
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
        warn("Reading from SBML model object is depricated, please do not use anymore")
        metabolites = parseMetaboliteInfoFromSBMLMod(model)
        _id = model.getId()
        return MeMoModel(metabolites=metabolites, _id=_id)

    def annotate(self) -> AnnotationResult:
        """Goes through the different bulk annotation methods and tries to annotate InChI strings to the metabolites
        in the model"""
        print(self._id)
        
        #TEMP FIX
        annotationDict = {"VMH" : annotateBiGG_id,
                          "ModelSEED": annotateModelSEED_id,
                          "BiGG": annotateBiGG_id
                          }

        origin_dbs = origin_databases(self.metabolites)
        origin_db = max(origin_dbs, key = lambda k: origin_dbs[k])
        print(f"ORIG DB {origin_db}")
        
        logger.debug(origin_db)
         # count the number of newly annotated metabolites
        anno_result = annotationDict[origin_db](self.metabolites)
        while True:
          old_res = AnnotationResult.fromAnnotation(anno_result)
          # BiGG
          temp_result = annotateBiGG(self.metabolites)
          #print("BiGG:",temp_result)
          anno_result = anno_result + temp_result
          # Use ChEBI
          temp_result = annotateChEBI(self.metabolites)
          #print("ChEBI:",temp_result)
          anno_result = anno_result + temp_result
          # GO BULK WISE ThORUGH BIGG AND VMH AND MODELSEED, try to extract as much as possible
          temp_result = annotateModelSEED(self.metabolites)
          #print("ModelSEED:", temp_result)
          anno_result = anno_result + temp_result
          #print("Total:", anno_result, "\n")
          if anno_result == old_res:
            break
          
        return anno_result


    def match(self, model2: MeMoModel, keep1ToMany:bool = True) -> pd.DataFrame:
        """ compares the metabolites of two models and returns a data frame with additional information """
        res_inchi = self.matchOnInchi(model2)
        res_db = self.matchOnDB(model2)
        res_name = self.matchOnName(model2)
        res = res_inchi.merge(res_db, how = "outer", on = ["met_id1","met_id2"],suffixes=["_inchi","_db"])
        res = res.merge(res_name, how = "outer", on = ["met_id1","met_id2"],suffixes=["","_name"])
        # TODO add comparison on the base of names
        return(res)

    
    def matchOnInchi(self, model2: MeMoModel) -> pd.DataFrame:
        # start with the comparison of inchi strings
        mod1_inchis = pd.DataFrame({"met_id" : [x.id for x in [y for y in self.metabolites]],
                "inchis" : [x._inchi_string for x in [y for y in self.metabolites]]})
        mod2_inchis = pd.DataFrame({"met_id" : [x.id for x in [y for y in model2.metabolites]],
                "inchis" : [x._inchi_string for x in [y for y in model2.metabolites]]})
        # go through the inchis of the second model and find the corresponding inchi in self
        matches = {"met_id1" : [],
                "met_id2" : [],
                "inchi_score":[],
                "inchi_string":[],
                "charge_diff" : []}


        
        mod1_inchis['Mol'] = mod1_inchis['inchis'].apply(inchiToMol)
        mod2_inchis['Mol'] = mod2_inchis['inchis'].apply(inchiToMol)


        mod1_inchis['fingerprint'] = mod1_inchis['Mol'].apply(molToRDK)
        mod2_inchis['fingerprint'] = mod2_inchis['Mol'].apply(molToRDK)

        mod1_inchis['normalized_inchi'] = mod1_inchis['Mol'].apply(molToNormalizedInchi)
        mod2_inchis['normalized_inchi'] = mod2_inchis['Mol'].apply(molToNormalizedInchi)

        for i in range(len(mod1_inchis)):
            inchi1 = mod1_inchis.loc[i, "inchis"]
            mol1   = mod1_inchis.loc[i, "Mol"]
            fp1 = mod1_inchis.loc[i, "fingerprint"]
            nminchi1 = mod1_inchis.loc[i, "normalized_inchi"]
            if inchi1 != None:
                id1 = mod1_inchis.loc[i,"met_id"]
                for j in range(len(mod2_inchis)):
                    inchi2 = mod2_inchis.loc[j, "inchis"]
                    mol2   = mod2_inchis.loc[i, "Mol"]
                    fp2 = mod2_inchis.loc[i, "fingerprint"]
                    nminchi2 = mod1_inchis.loc[i, "normalized_inchi"]
                    if inchi2 != None:
                        id2 = mod2_inchis.loc[j,"met_id"]
                        res = matchMetsByInchi(nminchi1, nminchi2, mol1, mol2, fp1, fp2)
                        if res[0] == True:
                            matches["met_id1"].append(id1)
                            matches["met_id2"].append(id2)
                            matches["inchi_string"].append(inchi1)
                            matches["charge_diff"].append(res[1])
                            matches["inchi_score"].append(1)
        inchiRes = pd.DataFrame(matches)
        return(inchiRes)

    def matchOnDB(self, model2: MeMoModel, threshold = 0, keep1ToMany = False) -> pd.DataFrame:
        # compare two models by the entries in the databases
        mets1 = self.metabolites
        mets2 = model2.metabolites
        results = {"met_id1":[],
                "met_id2":[],
                "DB_score":[],
                "charge_diff":[],
                "inchi_string":[]}
        for met1 in mets1:
            for met2 in mets2:
                jaccard = matchMetsByDB(met1,met2)
                if jaccard > threshold:
                    results["met_id1"].append(met1.id)
                    results["met_id2"].append(met2.id)
                    results["DB_score"].append(jaccard)
                    # try to calculate charge differences and add them as information
                    try:
                        charge_diff = met1._charge - met2._charge
                    except:
                        charge_diff = None
                    results["charge_diff"].append(charge_diff)
                    # try to find an inchi string for the pair and keep only the reference inchi
                    if met1._inchi_string != None or met2._inchi_string != None:
                        if met1._inchi_string != None:
                            inchi = met1._inchi_string
                        else:
                            inchi = met2._inchi_string
                        results["inchi_string"].append(inchi)
                    else:
                        results["inchi_string"].append(None)


        # remove one to many matches
        results = pd.DataFrame(results)
        if keep1ToMany== False : 
            results = results.sort_values(by = "DB_score", ascending = False)
            results = results.drop_duplicates("met_id1")
            results = results.sort_values(by = "DB_score", ascending = False)
            results = results.drop_duplicates("met_id2")
        return(results)

    def matchOnName(self, model2: MeMoModel, threshold = 0.6, keep1ToMany = False) -> pd.DataFrame:
        # compare two models by the entries in the databases
        mets1 = self.metabolites
        mets2 = model2.metabolites
        results = {"met_id1":[],
                "met_id2":[],
                "Name_score":[],
                "charge_diff":[],
                "inchi_string":[]}
        for met1 in mets1:
            for met2 in mets2:
                levenshtein = matchMetsByName(met1,met2)
                if levenshtein > threshold:
                    results["met_id1"].append(met1.id)
                    results["met_id2"].append(met2.id)
                    results["Name_score"].append(levenshtein)
                    # try to calculate charge differences and add them as information
                    try:
                        charge_diff = met1._charge - met2._charge
                    except:
                        charge_diff = None
                    results["charge_diff"].append(charge_diff)
                    # try to find an inchi string for the pair and keep only the reference inchi
                    if met1._inchi_string != None or met2._inchi_string != None:
                        if met1._inchi_string != None:
                            inchi = met1._inchi_string
                        else:
                            inchi = met2._inchi_string
                        results["inchi_string"].append(inchi)
                    else:
                        results["inchi_string"].append(None)


        # remove one to many matches
        results = pd.DataFrame(results)
        if keep1ToMany== False : 
            results = results.sort_values(by = "Name_score", ascending = False)
            results = results.drop_duplicates("met_id1")
            results = results.sort_values(by = "Name_score", ascending = False)
            results = results.drop_duplicates("met_id2")
        return(results)




