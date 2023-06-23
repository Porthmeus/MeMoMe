#Class:Metabolite
#values:
#    - Attributes {...}
#methods:
#    - annotate(self): try to get the InChI with different methods
#    - compare(self, Metabolite): -> [bool, score, method] compare itself to other metabolite classes
#         - compare by InChI
#         - compare by DB _ids
#         - compare by fuzzy name matching
#    - to_dict(self): forms a dictionary with all the information encoded in the metabolite
#    - get_id(self): returns the id
#    - get_inchi_string(self): returns the inchi_string
#    - get_annotations(self): returns a dictionary with the metabolite annotations
#
#TODO: add a function that forms a dataframe from the metabolite

from __future__ import annotations
from typing import List
from json import load
from numpy import NaN as npNaN
import warnings
from src.nameHandling import removeMetabolitePrefixSuffix
            
class MeMoMetabolite():
            
    """Class for "original" as well as "inferred"  information about a metabolite.
    The original information comes from a model, while the "inferred" information is either obtained though elaboration of the
    available information such as Inchi string or retrieved from the databases
    Metabolite is a class for holding directly know and inferred information regarding a metabolite 

    Parameters
    ----------
    _id : str
        the _id of the MeMoMetabolite - usually the core metabolite _id from the model.
    orig_ids : set()
        the original ID as specified in the model where the metabolite comes from, including the "M_" prefix and the compartment suffix:
    _model_id : str
        the ID of the model where the metabolite comes from
    names : set() 
        A list of common names and synonyms of the metabolite, will be converted to a set, to remove duplicates
    _inchi_string : string
       The Inchi string of the metabolite
    _formula : str
        The _formula of the metabolite
    _charge : float
       The _charge number of the metabolite
    annotations: dict()
       A dictionary whose keys are the names of the annotations, and the items are the database-specific metabolite _ids
    """

    def __init__(
        self,
        _id: str = None,
        orig_ids:list[str] = None,
        _model_id:str = None,
        names:list[str] = [],
        _inchi_string:str = None,
        _formula:str = None,
        _charge:int = None,
        annotations:dict = None,
    ) -> None:
        """Initialize MeMoMetaboblite

        Parameters
        ----------
        _id : str
            the _id of the MeMoMetabolite - usually the core metabolite _id from the model.
        orig_ids : set()
            the original ID as specified in the model where the metabolite comes from, including the "M_" prefix and the compartment suffix:
        _model_id : str
            the ID of the model where the metabolite comes from
        names : set() 
            A list of common names and synonyms of the metabolite, will be converted to a set, to remove duplicates
        _inchi_string : string
           The Inchi string of the metabolite
        _formula : str
            The _formula of the metabolite
        _charge : float
           The _charge number of the metabolite
        annotations: dict()
           A dictionary whose keys are the names of the annotations, and the items are the database-specific metabolite _ids
        """
        if _id != None:
            self._id =removeMetabolitePrefixSuffix(_id)
        else:
            self._id = _id
        self.orig_ids = set(orig_ids)
        self._model_id = _model_id
        self.names = set(names)
        self._formula = _formula
        self._inchi_string = _inchi_string
        self._charge = _charge
        self.annotations = annotations
    
    def set_id(self, new_id: str) -> None:
        ''' 
        set function for _id
        removes metabolite suffix and prefix
        '''
        if self._id != None:
            warnings.warn("changed metbolite _id from {old} to {new}".format(old = self._id, new = new_id))
        self._id = removeMetabolitePrefixSuffix(new_id)

    def set_model_id(self, new_model_id: str) -> None:
        ''' set function for _model_id '''
        if self._model_id != None:
            warnings.warn("changed metbolite _model_id from {old} to {new}".format(old = self._model_id, new = new_model_id))
        self._model_id = new_model_id
    
    def set_names(self, new_names: list[str]) -> None:
        ''' set function for names '''
        if self.names != []:
            warnings.warn("changed metbolite names from {old} to {new}".format(old = str(self.names), new = str(new_names)))
        self.names = set(new_names)
    
    def add_names(self, new_names: list[str]) -> None:
        ''' append new names to the list of metabolite names'''
        for x in new_names:
            self.names.add(x)

    def set_orig_ids(self, new_orig_ids: list[str]) -> None:
        ''' set function for orig_ids '''
        if self.orig_ids != []:
            warnings.warn("changed metbolite orig_ids from {old} to {new}".format(old = str(self.orig_ids), new = str(new_orig_ids)))
        self.orig_ids = set(new_orig_ids)
    
    def add_orig_ids(self, new_orig_ids: list[str]) -> None:
        ''' append new orig_ids to the list of metabolite orig_ids'''
        for x in new_orig_ids:
            self.orig_ids.add(x)
       
    def set_formula(self, new_formula: str) -> None:
        ''' set function for _formula '''
        if self._formula != None:
            warnings.warn("changed metbolite _formula from {old} to {new}".format(old = self._formula, new = new_formula))
        self._formula = new_formula

    def set_inchi_string(self, new_inchi_string: str) -> None:
        ''' set function for _inchi_string '''
        if self._inchi_string != None:
            warnings.warn("changed metbolite _inchi_string from {old} to {new}".format(old = self._inchi_string, new = new_inchi_string))
        self._inchi_string = new_inchi_string

    def set_charge(self, new_charge: int) -> None:
        ''' set function for _charge '''
        if self._charge != None:
            warnings.warn("changed metbolite _charge from {old} to {new}".format(old = str(self._charge), new = str(new_charge)))
        self._charge = new_charge

    def set_annotations(self, new_annotations: dict) -> None:
        ''' set function for annotations '''
        if self.annotations != None:
            warnings.warn("changed metbolite annotations from {old} to {new}".format(old = str(self.annotations), new = str(new_annotations)))
        self.annotations = new_annotations
    
    def add_annotations(self, new_annotations: dict) -> None:
        ''' append new annotations to the dict of metabolite annotations'''
        for x in new_annotations.keys():
            if x in self.annotations.keys():
                self.annotations[x].extend(new_annotations[x])
            else:
                self.annotations[x] = new_annotations[x]

    def get_unique_attributes(self) -> list:
        ''' Get a list of unique model attributes '''
        return([
            self._id,
            self._model_id,
            self._formula,
            self._inchi_string,
            self._charge
            ])

    def get_id(self) -> str:
        ''' Get the id '''
        return self._id

    def get_inchi_string(self) -> str:
        ''' Get the inchi string '''
        return self._inchi_string

    def get_annotations(self) -> dict:
        ''' Get the annotations '''
        return self.annotations

    def to_dict(self) -> dict:
        met_dict = dict()
        met_dict["id"] = self._id
        met_dict["orig_ids"] = self.orig_ids.copy()
        met_dict["model_id"] = self._model_id
        met_dict["names"] = self.names.copy()
        met_dict["inchi_string"] = self._inchi_string
        met_dict["formula"] = self._formula
        met_dict["charge"] = self._charge
        met_dict["annotations"] = self.annotations.copy() 
        return met_dict

    def merge(self, new_metabolite:MeMoMetabolite, keep_entries:bool = True, force:bool = False) -> None:
        ''' Merge two MeMoMetabolites into one entry 
        variables:
            new_metabolite:MeMoMetabolite - a MeMoMetabolite to merge into the self metabolite
            keep_entries:bool - if the unique attributes of self and the new_metabolite differ, which attribute should be used? Keep the original ones (default, True) or take the _id from the new_metabolite (False). Only used if force = True.
            force:bool - if unique attributes differ, should the merging be forced (True). Default: False, will raise an Error if attributes differ.
        '''
        
        # check if the IDs are the same, if not check what to do:
        if all([x == y for x,y in zip(self.get_unique_attributes(), new_metabolite.get_unique_attributes())]) != True:
            print([x == y for x,y in zip(self.get_unique_attributes(), new_metabolite.get_unique_attributes())])
            if force == False:
                raise ValueError("Unique attributes differ in the two metabolites {old} vs. {new}".format(old = str(self.get_unique_attributes()), new = str(new_metabolite.get_unique_attributes())))
            else:
                warnings.warn("Will merge, although unique attributes differ in the two metabolites {old} vs. {new}".format(old = str(self.get_unique_attributes()), new = str(new_metabolite.get_unique_attributes())))

        # if we want to overwrite the unique attributes do it here
        if keep_entries == False:
            self.set_id(new_metabolite._id)
            self.set_model_id(new_metabolite._model_id)
            self.set_formula(new_metabolite._formula)
            self.set_inchi_string(new_metabolite._inchi_string)
            self.set_charge(new_metabolite._charge)

        # add the non-unique attrbutes
        self.add_annotations(new_metabolite.annotations)
        self.add_names(list(new_metabolite.names))
        self.add_orig_ids(list(new_metabolite.orig_ids))


    def annotate(self):
        if self._inchi_string is None:
            #TODO: for now I am assuming that the download already happened, we can either add it on some external class or perfrom it here
            #downl_status = download()
            #if downl_status:
            #    print("Succesfully downloaded")

            # series of if/else statements to handle querying all the types of databases we have to retrieve the inchi string
            #we could define a sequence of databases/annotation from most relevant to less relevant and proceed following that sequence
            if "vmh" in self.annotations.keys():
                vmh_file = open("Databases/vmh.json","r")
                vmh_db = load(vmh_file)
                vmh_file.close()
                # access the vmh_db query results and saves the indexes associated to each _id
                abbreviations_index = {vmh_db["results"][i]["abbreviation"]:i for i in range(len(vmh_db["results"]))}
                vmh_id = self.annotations["vmh"]
                met_vmh_index = abbreviations_index[vmh_id]
                #TODO: add exception handling
                inchiString = vmh_db["results"][met_vmh_index]["inchiString"]
                if inchiString is not None and inchiString != npNaN and inchiString != "":
                    self._inchi_string = inchiString
            elif "bigg.metabolite" in self.annotations.keys():
                met_db_id = self.annotations["bigg.metabolite"]
                bigg_db = pd.read_csv("../Databases/BiGG.tsv", sep='\t',index_col=0) 
                #TODO: continue BiGG handling and other databases...
            
                




        # series of if/else statements to handle all the types of annotations we have (can be handled by calling the functions of conversion between i.e. smiles to inchi string
            

        return None
    
    def compare(self, metabolite):

        return None
