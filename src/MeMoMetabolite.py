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
from copy import deepcopy

from src.handle_metabolites_prefix_suffix import handle_metabolites_prefix_suffix, validateInchi
from src.annotation.annotateInchiRoutines import findOptimalInchi

class MeMoMetabolite():
    """Class for "original" as well as "inferred"  information about a metabolite.
    The original information comes from a model, while the "inferred" information is either obtained though elaboration of the
    available information such as Inchi string or retrieved from the databases
    Metabolite is a class for holding directly know and inferred information regarding a metabolite 
    """

    def __init__(
            self,
            _id: str|None = None,
            orig_ids: list[str]|None = None,
            _model_id: str|None = None,
            names: list[str]|None = None,
            names_source: list[str]|None = None,
            _inchi_string: str|None = None,
            _inchi_source: str|None = None,
            _formula: str|None = None,
            _charge: int|None = None,
            # TODO _charge_source: str|None = None,
            _pKa: dict[str,float]|None= None,
            # TODO  _pKa_source: dict|None = None,
            _pKb: dict[str,float]|None = None,
            # TODO _pKb_source: dict|None = None,
            annotations: dict[str,list[str]]|None = None,
            annotations_source: dict[str,list[str]]|str|None = None,
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
            self._id = None 
    
        
        # use setter function to acutally set the values to make sure to remove
        # duplicates and sort all the lists - this is important and you should
        # not circumvent this step just for fun
        if orig_ids is None:
            self.orig_ids = []
        else:
            self.orig_ids = []
            self.set_orig_ids(orig_ids)

        if _model_id is None:
            self._model_id = None
        else:
            self._model_id = None
            self.set_model_id(_model_id)

        if names is None:
            self.names = []
            self.names_source =[]
        else:
            self.names = []
            self.names_source =[]
            if names_source is None:
                names_source = "model"
            self.set_names(names, source = names_source)

        if _formula is None:
            self._formula = None
        else:
            self._formula = None
            self.set_formula(_formula)

        if _inchi_string is None:
            self._inchi_string = None
            self._inchi_source = None
        else:
            self._inchi_string = None
            if _inchi_source is None:
                _inchi_source = "model"
            self.set_inchi_string(_inchi_string, source = _inchi_source)

        if _charge is None:
            self._charge = None
        else:
            self._charge = None
            self.set_charge(_charge)

        if _pKa is None:
            self._pKa = {}
        else:
            self._pKa = {}
            self.set_pKa(_pKa)
        
        if _pKb is None:
            self._pKb = {}
        else:
            self._pKb = {}
            self.set_pKb(_pKb)

        if annotations is None:
            self.annotations = {}
            self.annotations_source = {}
        else:
            self.annotations = {}
            self.annotations_source = {}
            if annotations_source is None:
                annotations_source = "model"
            self.set_annotations(annotations, source = annotations_source)


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

    def set_names(self, new_names: list[str], source:str|list[str]|None=None) -> None:
        """ set function for names """
        if source == None:
            source = self.names_source
        elif type(source) == str:
            source = [source]*len(new_names)

        if len(source) != len(new_names):
            raise ValueError("List of new names have a different length than the one for the sources. Consider providing a corresponding soure list to the new names via the source parameter.")

        if self.names:
            warnings.warn("changed metbolite names from {old} to {new}".format(old=str(self.names), new=str(new_names)))
        
        # we want the names to be sorted lexographically, thus we have to sort the sources accordingly, so that the index of the name corresponds to the index of the source. This is done with this rather complicated line below (stackoverflow for the rescue)
        new_names, source = zip(*sorted(zip(new_names,source))) # sort both 
        self.names = list(new_names)
        self.names_source = list(source)

    def add_names(self, new_names: list[str], source: str) -> int:
        ''' append new names to the list of metabolite names'''
        old_names = deepcopy(self.names)
        self.names.extend(new_names)
        # remove duplicates and sort lexographically
        self.names = list(dict.fromkeys(self.names)) # use dict instead of set to preserve order
        # we allow addition of names only by one source at a time - however, if we want to add several names we need to extend the source list by the same number of elements as the name list
        self.names_source.extend([source]*(len(self.names)-len(old_names))) # extend source list
        # we want the names to be sorted lexographically, thus we have to sort the sources accordingly, so that the index of the name corresponds to the index of the source. This is done with this rather complicated line below (stackoverflow for the rescue)
        new_names, new_sources = zip(*sorted(zip(self.names,self.names_source))) # sort both list for self.names
        self.names = list(new_names)
        self.names_source = list(new_sources)
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
        old_orig_ids = deepcopy(self.orig_ids)
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

    def set_inchi_string(self, new_inchi_string: str, source: str) -> int:
        """ set function for _inchi_string """
        old_inchi = deepcopy(self._inchi_string)
        
        # validate inchi string - handles also empty strings
        new_inchi_string = validateInchi(new_inchi_string)
        
        if self._inchi_string is not None:
            warnings.warn("changed metbolite _inchi_string from {old} to {new}".format(old=self._inchi_string,
                                                                                       new=new_inchi_string))
        self._inchi_string = new_inchi_string
        self._inchi_source = source
        return int(self._inchi_string != old_inchi)
    
    def add_inchi_string(self, new_inchi_string:str, source:str) -> int:
        """ compares to inchis and takes the most appropiate one"""
        changed = 0
        
        # handle invalid inchi strings -handles also empty strings
        new_inchi_string = validateInchi(new_inchi_string)
        if new_inchi_string == None:
            return changed
        
        if self._inchi_string != None:
            old = self._inchi_string
            new = findOptimalInchi([self._inchi_string, new_inchi_string], charge = self._charge)
            if new is None:
                raise NotImplementedError()
            if new != old:
                #print(self.id) # nice feature for debugging
                #print("OLD:",old)
                #print("NEW:",new)
                self._inchi_string = new
                self._inchi_source = source
                changed =  1
        else:
            self._inchi_string = new_inchi_string
            self._inchi_source = source
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
        old_pKa = deepcopy(self._pKa)
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
        old_pKb = deepcopy(self._pKb)
        if self._pKb is not None and self._pKb != {}:
            for key,value in new_pKb.items():
                if key in self._pKb.keys():
                    self._pKb[key] = (self._pKb[key] + value)/2
                else:
                    self._pKb[key] = value
        return int(old_pKb != self._pKb)

    def set_pKb(self, new_pKb:dict) -> None:
        """ set function for _pKb """
        if self._pKb is not None and self._pKb != {}:
            old = ";".join([x+":"+str(y) for x,y in self._pKb.items()])
            new = ";".join([x+":"+str(y) for x,y in new_pKb.items()])
            warnings.warn("changed metabolite _pKs from {old} to {new}".format(old = old, new =new))
        self._pKb = new_pKb
    
    def set_annotations(self, new_annotations: dict, source:str|dict|None=None) -> None:
        """ set function for annotations """

        # try to make sources if there is not much information
        if source == None:
            source == deepcopy(self.annotations_source)
        elif type(source) == str:
            val = source
            source = dict()
            for x,y in new_annotations.items():
                source[x] = [val]*len(y)

        # first check if the structure of annotations and sources are the same
        if new_annotations.keys() != source.keys():
            raise ValueError('''Keys of the new annotation dictionary differs from the keys in source. 
            Annotation keys: {key1}
            Source keys: {key2}'''.format(key1 = new_annotations.keys(), key2 = source.keys()))
        elif [len(x) for x in new_annotations.values()] != [len(y) for y in source.values()]:
            raise ValueError('''Entries length differs for new annotation and source dictionary. Consider to provide a source dictionary via the source value''')

        # check if there is an old annotation which will be overwritten
        if self.annotations != {} and self.annotations != None:
            warntrigger = True
            old_annotations = deepcopy(self.annotations)
        else:
            warntrigger = False
        
        # delete the old information
        self.annotations = dict()
        self.annotations_source = dict()
        
        # sort both new dictionaries and remove duplicates 
        for x in new_annotations.keys():
            new_annolst = []
            new_sourcelst = []
            for anno,src in zip(new_annotations[x], source[x]):
                if not anno in new_annolst:
                    new_annolst.append(anno)
                    new_sourcelst.append(src)
            new_annolst,new_sourcelst = zip(*sorted(zip(new_annolst,new_sourcelst)))
            self.annotations[x] = list(new_annolst)
            self.annotations_source[x] = list(new_sourcelst)
        
        # trigger the warning and give some information to the user
        if warntrigger == True:
            warnings.warn("changed metbolite annotations from {old} to {new}".format(
                new=str(self.annotations),
                old=str(old_annotations)))

    def add_annotations(self, new_annotations: dict, source:str) -> int:
        """ append new annotations to the dict of metabolite annotations
        check if there have been actually added new annotation, if so return 1,else 0"""
        old_annotation = deepcopy(self.annotations)
        for x in new_annotations.keys():
            # check if there is already annotations with that key
            if x in self.annotations.keys():
                self.annotations[x].extend(new_annotations[x])
                # remove duplicates and sort list
                new_annolst = list(dict.fromkeys(self.annotations[x])) # hacky way to create uniques
                new_sources = self.annotations_source[x]
                # here we extend the list number of source by the number of new annotations
                new_sources.extend([source]*(len(new_annolst) - len(old_annotation[x])))
                # we want the annotations to be sorted lexographically, thus we have
                # to sort the sources accordingly, so that the index of the
                # anno corresponds to the index of the source. This is done
                # with this rather complicated line below (stackoverflow for
                # the rescue)
                new_annolst, new_sources = zip(*sorted(zip(new_annolst,new_sources)))
                self.annotations[x] = list(new_annolst)
                self.annotations_source[x] = list(new_sources)
            else: # otherwise create a new dictionary entry
                # remove duplicates and sort list
                new_annolst = list(dict.fromkeys(new_annotations[x])) # hacky way to create uniques
                # here we extend the list number of source by the number of new annotations
                new_sources = [source]*len(new_annolst)
                # we want the annotations to be sorted lexographically, thus we have
                # to sort the sources accordingly, so that the index of the
                # anno corresponds to the index of the source. This is done
                # with this rather complicated line below (stackoverflow for
                # the rescue)
                new_annolst,new_sources = zip(*sorted(zip(new_annolst,new_sources)))
                self.annotations[x] = list(new_annolst)
                self.annotations_source[x] = list(new_sources)
        
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
            self.set_inchi_string(new_metabolite._inchi_string, new_metabolite._inchi_source)
            self.set_charge(new_metabolite._charge)

        # add the non-unique attrbutes including their sources of annotation
        for key in new_metabolite.annotations.keys():
            for i in range(len(new_metabolite.annotations[key])):
                new_anno_dic = {key : [new_metabolite.annotations[key][i]]}
                new_source = new_metabolite.annotations_source[key][i]
                self.add_annotations(new_anno_dic, source = new_source)

        # need to iterate, because source is only accepted as single string as
        # variable of add_names(), yet it can contain several different values
        # if derived from another model
        for i in range(len(new_metabolite.names)):
            self.add_names([new_metabolite.names[i]], new_metabolite.names_source[i])
        self.add_orig_ids(list(new_metabolite.orig_ids))

    def annotate(self):
#        if self._inchi_string is None:
            # TODO: for now I am assuming that the download already happened,
            # we can either add it on some external class or perfrom it here
            # downl_status = download()
            # if downl_status:
            #    print("Succesfully downloaded")

            # series of if/else statements to handle querying all the types of
            # databases we have to retrieve the inchi string
            # we could define a sequence of databases/annotation from most
            # relevant to less relevant and proceed following that sequence
        warnings.warn("MeMoMetabolite.annotate() is not implemented yet")
        return None

    def compare(self, metabolite):

        warnings.warn("MeMoMetabolite.compare() is not implemented yet")
        return None

    @property
    def id(self):
        return self._id
