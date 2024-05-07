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
import logging
from src.annotateChEBI import annotateChEBI
from src.annotateBiGG import annotateBiGG, annotateBiGG_id
from src.annotateModelSEED import annotateModelSEED, annotateModelSEED_id
from src.annotateVMH import annotateVMH, annotateVMH_id
from src.annotateAux import AnnotationResult
from src.matchMets import matchMetsByDB, matchMetsByInchi, matchMetsByName, NeutraliseCharges
from src.parseMetaboliteInfos import parseMetaboliteInfoFromSBML, parseMetaboliteInfoFromSBMLMod, \
    parseMetaboliteInfoFromCobra
from src.MeMoMetabolite import MeMoMetabolite
from src.annotateInchiRoutines import inchiToMol, molToRDK, molToNormalizedInchi
from src.origin_databases import origin_databases
from rdkit import Chem
import logging

logger = logging.getLogger('logger')


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
    def fromPath(cls, sbmlfile: Path) -> MeMoModel :
        """ Read the model from the SBML file """
        metabolites = parseMetaboliteInfoFromSBML(sbmlfile, validate=True)
        cobra_model, errors = cb.io.sbml.validate_sbml_model(sbmlfile)
        if any([len(x) > 0 for x in errors.values() ]):
          logger.error(f"There were problems with the sbml model {sbmlfile}")
          logger.error(f"{errors}")
        _id = cobra_model.id
        return MeMoModel(metabolites=metabolites, cobra_model=cobra_model, _id=_id)

    @classmethod
    def fromModel(cls, model: cb.Model) -> MeMoModel:
        """ Read the model from the cobra model """
        cobra_model = model
        _id = model.id
        metabolites = parseMetaboliteInfoFromCobra(model)
        return MeMoModel(metabolites=metabolites,cobra_model=cobra_model, _id=_id)

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
        annotationDict = {"VMH" : annotateVMH_id,
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
          print("BiGG:",temp_result)
          anno_result = anno_result + temp_result
          # Use ChEBI
          temp_result = annotateChEBI(self.metabolites)
          print("ChEBI:",temp_result)
          anno_result = anno_result + temp_result
          temp_result = annotateModelSEED(self.metabolites)
          anno_result = anno_result + temp_result
          print("ModelSEED:", temp_result)
          temp_result = annotateVMH(self.metabolites)
          print("VMH:", temp_result)
          anno_result = anno_result + temp_result
          print("Total:", anno_result, "\n")
          print("\n\n\n\n\n\n\n\n\n")
          if anno_result == old_res:
            break
          
        return anno_result



    def match(self, model2: MeMoModel, keep1ToMany:bool = True, output_names: bool = False, output_dbs: bool = False, keepUnmatched: bool = False) -> pd.DataFrame:
        """ compares the metabolites of two models and returns a data frame with additional information 
        output_names: If true, output the names of the metabolites that led to match based on levenshtein
        """
        res_inchi = self.matchOnInchi(model2, keep1ToMany = keep1ToMany)
        res_db = self.matchOnDB(model2, output_dbs = output_dbs, keep1ToMany= keep1ToMany)
        res_name = self.matchOnName(model2, output_names = output_names, keep1ToMany= keep1ToMany)

        res = res_inchi.merge(res_db, how = "outer", on = ["met_id1","met_id2"],suffixes=["_inchi","_db"])
        res = res.merge(res_name, how = "outer", on = ["met_id1","met_id2"],suffixes=["","_name"])

        # add the unmatched metabolites
        if keepUnmatched == True:
            miss_mets1 = list(set([x.id for x in self.metabolites]) - set(res["met_id1"]))
            missing_df1 = pd.DataFrame({"met_id1":miss_mets1})
            miss_mets2 = list(set([x.id for x in model2.metabolites]) - set(res["met_id2"]))
            missing_df2 = pd.DataFrame({"met_id2":miss_mets2})
            res = pd.concat([res,missing_df1, missing_df2])

        return(res)

    
    def matchOnInchi(self, model2: MeMoModel, keep1ToMany:bool = False) -> pd.DataFrame:

        # create data frames containing the information for comparison 
        mod1_inchis = pd.DataFrame({"met_id" : [x.id for x in [y for y in self.metabolites]],
                "inchis" : [x._inchi_string for x in [y for y in self.metabolites]]})
        mod2_inchis = pd.DataFrame({"met_id" : [x.id for x in [y for y in model2.metabolites]],
                "inchis" : [x._inchi_string for x in [y for y in model2.metabolites]]})

        # do some precalculations for speed up
        # precalculate the mol representation for the inchi in rdkit
        mod1_inchis['Mol'] = mod1_inchis['inchis'].apply(inchiToMol)
        mod2_inchis['Mol'] = mod2_inchis['inchis'].apply(inchiToMol)
        
        # precalculate the fingerprints
        mod1_inchis['fingerprint'] = mod1_inchis['Mol'].apply(molToRDK)
        mod2_inchis['fingerprint'] = mod2_inchis['Mol'].apply(molToRDK)
        
        # normalize the inchi strings
        mod1_inchis['normalized_inchi'] = mod1_inchis['Mol'].apply(molToNormalizedInchi)
        mod2_inchis['normalized_inchi'] = mod2_inchis['Mol'].apply(molToNormalizedInchi)

        # remove charges
        mask1 = ~pd.isna(mod1_inchis.Mol)
        if sum(mask1) > 0:
            mod1_inchis.loc[mask1,'neutralized_charge_mol'] = mod1_inchis.loc[mask1,'Mol'].apply(NeutraliseCharges)
        else:
            mod1_inchis["neutralized_charge_mol"] = None

        mask2 = ~pd.isna(mod2_inchis.Mol)
        if sum(mask2) > 0:
            mod2_inchis.loc[mask2,'neutralized_charge_mol'] = mod2_inchis.loc[mask2,'Mol'].apply(NeutraliseCharges)
        else:
            mod2_inchis["neutralized_charge_mol"] = None

        # go through the inchis of the second model and find the corresponding inchi in self
        matches = {"met_id1" : [],
                "met_id2" : [],
                "inchi_score":[],
                "inchi_string":[],
                "charge_diff" : []}

        # loop over the first models metabolites
        for i in range(len(mod1_inchis)):
            # get the first inchi
            inchi1 = mod1_inchis.loc[i, "inchis"]
            # check if there is acutally an inchi, or whether the metabolite has none
            if inchi1 != None:
                # assign the precalculated values for the inchi
                mol1   = mod1_inchis.loc[i, "Mol"]
                fp1 = mod1_inchis.loc[i, "fingerprint"]
                nminchi1 = mod1_inchis.loc[i, "normalized_inchi"]
                ntchrmol1 = mod1_inchis.loc[i, "neutralized_charge_mol"]
                id1 = mod1_inchis.loc[i,"met_id"]

                # loop through the metabolites of the second model
                for j in range(len(mod2_inchis)):
                    inchi2 = mod2_inchis.loc[j, "inchis"]
                    # check if there is acutally an inchi, or whether the metabolite has none
                    if inchi2 != None:
                        # assign the precalculated values for the inchi
                        mol2   = mod2_inchis.loc[j, "Mol"]
                        fp2 = mod2_inchis.loc[j, "fingerprint"]
                        nminchi2 = mod2_inchis.loc[j, "normalized_inchi"]
                        ntchrmol2 = mod2_inchis.loc[j, "neutralized_charge_mol"]
                        id2 = mod2_inchis.loc[j,"met_id"]
                        
                        # do the actual matching
                        res = matchMetsByInchi(nminchi1, nminchi2, mol1, mol2, fp1, fp2, ntchrmol1, ntchrmol2)
                        
                        # write the results into the data frame
                        if res[0] == True or keep1ToMany == True:
                            matches["met_id1"].append(id1)
                            matches["met_id2"].append(id2)
                            matches["inchi_string"].append(inchi1)
                            matches["charge_diff"].append(res[1])
                            matches["inchi_score"].append(int(res[0]))
        inchiRes = pd.DataFrame(matches)
        return(inchiRes)

    def matchOnDB(self, model2: MeMoModel, threshold = 0, keep1ToMany = False, output_dbs: bool = False ) -> pd.DataFrame:
        print(f"Output_dbs set to {output_dbs}")
        # compare two models by the entries in the databases

        # get the metabolites and set up a results data frame
        mets1 = self.metabolites
        mets2 = model2.metabolites
        results = {"met_id1":[],
                "met_id2":[],
                "DB_score":[],
                "charge_diff":[],
                "inchi_string":[],
                "commonDBs" : [],
                "commonIds": [],
                "allIds": []}
        
        # if output_dbs should not be reported, remove the entries here from the results data frame
        if not output_dbs:
          results.pop("commonDBs")
          results.pop("commonIds")
          results.pop("allIds")

        
        # loop through all metabolites
        for met1 in mets1:
            for met2 in mets2:
                jaccard = matchMetsByDB(met1,met2)
                if jaccard.score > threshold:
                    results["met_id1"].append(met1.id)
                    results["met_id2"].append(met2.id)
                    results["DB_score"].append(jaccard.score)
                    if output_dbs:
                        results["commonDBs"] = jaccard.commonDBs
                        results["commonIds"] = jaccard.commonIds
                        results["allIds"] = jaccard.allIds
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

    def matchOnName(self, model2: MeMoModel, threshold = 0.6, keep1ToMany = False, output_names: bool = False) -> pd.DataFrame:
        print(f"Output_names set to {output_names}")
        # compare two models by the entries in the databases
        mets1 = self.metabolites
        mets2 = model2.metabolites
        results = {"met_id1":[],
                "met_id2":[],
                "Name_score":[],
                "charge_diff":[],
                "inchi_string":[],
                   "name_id1" : [],
                   "name_id2": []}
        

        if not output_names:
          results.pop("name_id1")
          results.pop("name_id2")
        for met1 in mets1:
            for met2 in mets2:
                levenshtein = matchMetsByName(met1,met2)
                if levenshtein.score > threshold:
                    results["met_id1"].append(met1.id)
                    results["met_id2"].append(met2.id)
                    if output_names:
                      results["name_id1"].append(levenshtein.name_id1)
                      results["name_id2"].append(levenshtein.name_id2)
                    results["Name_score"].append(levenshtein.score)
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




