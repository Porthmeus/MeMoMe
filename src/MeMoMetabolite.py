# Class:Metabolite
# values:
#    - Attributes {...}
# methods:
#    - annotate(self): try to get the InChI with different methods
#    - compare(self, Metabolite): -> [bool, score, method] compare itself to other metabolite classes
#         - compare by InChI
#         - compare by DB _ids
#         - compare by fuzzy name matching

#### TODO:
# - add a function which forms a dataframe from the metabolite
# - getter and setter functions should ensure sorted and unique lists and dictionaries
# - add_ functions should return 1 if really something was added, else it should return 0

from __future__ import annotations

import warnings
from json import load

import pandas as pd
from numpy import NaN as npNaN

from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix


class MeMoMetabolite():
    """Class for "original" as well as "inferred"  information about a metabolite.
    The original information comes from a model, while the "inferred" information is either obtained though elaboration of the
    available information such as Inchi string or retrieved from the databases
    Metabolite is a class for holding directly know and inferred information regarding a metabolite 
    """

    def __init__(
            self,
            _id: str = None,
            orig_ids: list[str] = [],
            _model_id: str = None,
            names=None,
            _inchi_string: str = None,
            _formula: str = None,
            _charge: int = None,
            _pKa: dict = {},
            _pKb: dict = {},
            annotations: dict = {},
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
        _pKa : float
            The pKa of the metabolite
        _pKb : float
            The pKb of the metabolite
        annotations: dict()
           A dictionary whose keys are the names of the annotations, and the items are the database-specific metabolite _ids
        """
        if _id is not None:
            self._id = handle_metabolites_prefix_suffix(_id)
        else:
            # TODO Why assign None? What else?
            self._id = None 

        self.orig_ids = []
        self._model_id = None
        self.names = [] 
        self._formula = None
        self._inchi_string = None
        self._charge = None
        self._pKa = {}
        self._pKb = {}
        self.annotations = {}
        
        # use setter function to acutally set the values to make sure to remove
        # duplicates and sort all the lists - this is important and you should
        # not circumvent this step just for fun
        self.set_orig_ids(orig_ids)
        self.set_model_id(_model_id)
        if names is None:
            names = []
        else:
            self.set_names(names)
        self.set_formula(_formula)
        self.set_inchi_string(_inchi_string)
        self.set_charge(_charge)
        self.set_pKa(_pKa)
        self.set_pKb(_pKb)
        self.set_annotations(annotations)


    def set_id(self, new_id: str) -> None:
        """
        set function for _id
        removes metabolite suffix and prefix
        """
        if self._id is not None:
            warnings.warn("changed metbolite _id from {old} to {new}".format(old=self._id, new=new_id))
        self._id = handle_metabolites_prefix_suffix(new_id)

    def set_model_id(self, new_model_id: str) -> None:
        """ set function for _model_id """
        if self._model_id is not None:
            warnings.warn(
                "changed metbolite _model_id from {old} to {new}".format(old=self._model_id, new=new_model_id))
        self._model_id = new_model_id

    def set_names(self, new_names: list[str]) -> None:
        """ set function for names """
        if self.names:
            warnings.warn("changed metbolite names from {old} to {new}".format(old=str(self.names), new=str(new_names)))
        new_names = list(set(new_names))
        new_names.sort()
        self.names = new_names

    def add_names(self, new_names: list[str]) -> int:
        ''' append new names to the list of metabolite names'''
        old_names = self.names.copy()
        for x in new_names:
            self.names.append(x)
        # remove duplicates and sort lexographically
        self.names = list(set(self.names))
        self.names.sort()
        # return if the names have changed
        return int(old_names != self.names)

    def set_orig_ids(self, new_orig_ids: list[str]) -> None:
        """ set function for orig_ids """
        if self.orig_ids != []:
            warnings.warn(
                "changed metbolite orig_ids from {old} to {new}".format(
                    old=str(self.orig_ids),
                    new=str(new_orig_ids)))
        self.orig_ids = list(set(new_orig_ids))
        self.orig_ids.sort()

    def add_orig_ids(self, new_orig_ids: list[str]) -> int:
        """ append new orig_ids to the list of metabolite orig_ids
        returns 1 if change has occured, else 0"""
        old_orig_ids = self.orig_ids.copy()
        for x in new_orig_ids:
            self.orig_ids.append(x)
        # make sure that there are no duplicates
        self.orig_ids = list(set(self.orig_ids))
        self.orig_ids.sort()
        # check whethter there was acutally an addition of values
        return int(old_orig_ids != self.orig_ids)

    def set_formula(self, new_formula: str) -> None:
        """ set function for _formula """
        if self._formula is not None:
            warnings.warn("changed metbolite _formula from {old} to {new}".format(old=self._formula, new=new_formula))
        self._formula = new_formula

    def set_inchi_string(self, new_inchi_string: str) -> int:
        """ set function for _inchi_string """
        old_inchi = self._inchi_string
        if self._inchi_string is not None:
            warnings.warn("changed metbolite _inchi_string from {old} to {new}".format(old=self._inchi_string,
                                                                                       new=new_inchi_string))
        self._inchi_string = new_inchi_string
        return int(self._inchi_string != old_inchi)
    
    def add_inchi_string(self, new_inchi_string:str) -> int:
        """ compares to inchis and takes the most appropiate one"""
        changed = 0
        if self._inchi_string != None:
            old = self._inchi_string
            new = findOptimalInchi(self._inchi_string, new_inchi_string)
            if new != old:
                self._inchi_string = new
                changed =  1
        else:
            self._inchi_string = new_inchi_string
            changed = 1
        return changed
    

    def set_charge(self, new_charge: int) -> None:
        """ set function for _charge """
        if self._charge is not None:
            warnings.warn(
                "changed metbolite _charge from {old} to {new}".format(old=str(self._charge), new=str(new_charge)))
        self._charge = new_charge
    
    def add_pKa(self, new_pKa:float) -> int:
        """ add a new pKa value, this is a simple mean calculation and might need some review"""
        # TODO: check for correct calulations of pKa mean
        old_pKa = self._pKa.copy()
        if self._pKa is not None:
            for key,value in new_pKa.items():
                if key in self._pKa.keys():
                    self._pKa[key] = (self._pKa[key] + value)/2
                else:
                    self._pKa[key] = value
        return int(old_pKa != self._pKa)


    def set_pKa(self, new_pKa:dict) -> None:
        """ set function for _pKa """
        if self._pKa is not None and self._pKa != {}:
            old = ";".join([x+":"+str(y) for x,y in self._pKa.items()])
            new = ";".join([x+":"+str(y) for x,y in new_pKa.items()])
            warnings.warn("changed metabolite _pKs from {old} to {new}".format(old = old, new =new))
        self._pKa = new_pKa

    def add_pKb(self, new_pKb:float) -> int:
        """ add a new pKb value, this is a simple mean calculation and might need some review"""
        # TODO: check for correct calulations of pKb mean
        old_pKb = self._pKb.copy()
        if self._pKb is not None and self._pKb != {}:
            for key,value in new_pKb.items():
                if key in self._pKb.keys():
                    self._pKb[key] = (self._pKb[key] + value)/2
                else:
                    self._pKb[key] = value
        return int(old_pKb != self._pKb)

    def set_pKb(self, new_pKb:dict) -> None:
        """ set function for _pKb """
        if self._pKb is not None:
            old = ";".join([x+":"+str(y) for x,y in self._pKb.items()])
            new = ";".join([x+":"+str(y) for x,y in new_pKb.items()])
            warnings.warn("changed metabolite _pKs from {old} to {new}".format(old = old, new =new))
        self._pKb = new_pKb
    
    def set_annotations(self, new_annotations: dict) -> None:
        """ set function for annotations """
        if self.annotations != {} and self.annotations != None:
            warnings.warn("changed metbolite annotations from {old} to {new}".format(
                old=str(self.annotations),
                new=str(new_annotations)))
        self.annotations = new_annotations

    def add_annotations(self, new_annotations: dict) -> int:
        """ append new annotations to the dict of metabolite annotations
        check if there have been actually added new annotation, if so return 1,else 0"""
        old_annotation = self.annotations.copy()
        for x in new_annotations.keys():
            if x in self.annotations.keys():
                self.annotations[x].extend(new_annotations[x])
                # remove duplicates and sort list
                new_annolst = list(set(self.annotations[x]))
                new_annolst.sort()
                self.annotations[x] = new_annolst
            else:
                # remove duplicates and sort list
                new_annolst = list(set(new_annotations[x]))
                new_annolst.sort()
                self.annotations[x] = new_annolst
        
        # check if the annotations have been changed
        return int(old_annotation!=self.annotations)

    def get_unique_attributes(self) -> list:
        """ Get a list of unique model attributes """
        return ([
            self._id,
            self._model_id,
            self._formula,
            self._inchi_string,
            self._charge
        ])

    def merge(self, new_metabolite: MeMoMetabolite, keep_entries: bool = True, force: bool = False) -> None:
        ''' Merge two MeMoMetabolites into one entry 
        variables:
            new_metabolite:MeMoMetabolite - a MeMoMetabolite to merge into the self metabolite
            keep_entries:bool - if the unique attributes of self and the new_metabolite differ, which attribute should be used? Keep the original ones (default, True) or take the _id from the new_metabolite (False). Only used if force = True.
            force:bool - if unique attributes differ, should the merging be forced (True). Default: False, will raise an Error if attributes differ.
        '''

        # check if the IDs are the same, if not check what to do:
        if all([x == y for x, y in zip(self.get_unique_attributes(), new_metabolite.get_unique_attributes())]) != True:
            print([x == y for x, y in zip(self.get_unique_attributes(), new_metabolite.get_unique_attributes())])
            if force == False:
                raise ValueError("Unique attributes differ in the two metabolites {old} vs. {new}".format(
                    old=str(self.get_unique_attributes()), new=str(new_metabolite.get_unique_attributes())))
            else:
                warnings.warn(
                    "Will merge, although unique attributes differ in the two metabolites {old} vs. {new}".format(
                        old=str(self.get_unique_attributes()), new=str(new_metabolite.get_unique_attributes())))

        # if we want to overwrite the unique attributes do it here
        if not keep_entries:
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
            # TODO: for now I am assuming that the download already happened,
            # we can either add it on some external class or perfrom it here
            # downl_status = download()
            # if downl_status:
            #    print("Succesfully downloaded")

            # series of if/else statements to handle querying all the types of
            # databases we have to retrieve the inchi string
            # we could define a sequence of databases/annotation from most
            # relevant to less relevant and proceed following that sequence
            if "vmh" in self.annotations.keys():
                vmh_file = open("Databases/vmh.json", "r")
                vmh_db = load(vmh_file)
                vmh_file.close()
                # access the vmh_db query results and saves the indexes associated to each _id
                abbreviations_index = {vmh_db["results"][i]["abbreviation"]: i for i in range(len(vmh_db["results"]))}
                vmh_id = self.annotations["vmh"]
                met_vmh_index = abbreviations_index[vmh_id]
                # TODO: add exception handling
                inchiString = vmh_db["results"][met_vmh_index]["inchiString"]
                if inchiString is not None and inchiString != npNaN and inchiString != "":
                    self._inchi_string = inchiString
            elif "bigg.metabolite" in self.annotations.keys():
                met_db_id = self.annotations["bigg.metabolite"]
                bigg_db = pd.read_csv("../Databases/BiGG.tsv", sep='\t', index_col=0)
                # TODO: continue BiGG handling and other databases...

        # series of if/else statements to handle all the types of annotations
        # we have (can be handled by calling the functions of conversion
        # between i.e. smiles to inchi string

        return None

    def compare(self, metabolite):

        return None

    @property
    def id(self):
        return self._id
