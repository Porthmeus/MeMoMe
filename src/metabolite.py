#Class:Metabolite
#values:
#    - Attributes {...}
#methods:
#    - annotate(self): try to get the InChI with different methods
#    - compare(self, Metabolite): -> [bool, score, method] compare itself to other metabolite classes
#         - compare by InChI
#         - compare by DB ids
#         - compare by fuzzy name matching

from typing import List
from json import load
from numpy import NaN as npNaN
            
class MeMoMetabolite():
            
    """Class for "original" as well as "inferred"  information about a metabolite.
    The original information comes from a model, while the "inferred" information is either obtained though elaboration of the
    available information such as Inchi string or retrieved from the databases
    Metabolite is a class for holding directly know and inferred information regarding a metabolite 

    Parameters
    ----------
    orig_id : str
        the original ID as specified in the model where the metabolite comes from
    model : str
        the ID of the model where the metabolite comes from
    names : list()
        A list of synonims of the metabolite
    inchi_string : string
       The Inchi string of the metabolite
    formula : str
        The formula of the metabolite
    charge : float
       The charge number of the metabolite
    annotations: dict()
       A dictionary whose keys are the names of the databases/types of annotations, and the items are the database(or annotations)-specific metabolite ids
    """

    def __init__(
        self,
        orig_id:str = None,
        model:str = None,
        names:List[str] = [],
        inchi_string:str = None,
        formula:str = None,
        charge:float = None,
        annotations:dict = None,
    ) -> None:
        """Initialize MeMoMetaboblite

        Parameters
        ----------
        orig_id : str
            the original ID as specified in the model where the metabolite comes from
        model : str
            the ID of the model where the metabolite comes from
        names : list()
            A list of synonims of the metabolite
        inchi_string : string
           The Inchi string of the metabolite
        formula : str
            The formula of the metabolite
        charge : float
           The charge number of the metabolite
        annotations: dict()
           A dictionary whose keys are the names of the annotations, and the items are the database-specific metabolite ids
        """
        self.orig_id = orig_id
        self.model = model
        self.names = names
        self.formula = formula
        self.inchi_string = inchi_string
        self.charge = charge
        self.annotations = annotations

        
    def annotate(self):
        if self.inchi_string is None:
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
                # access the vmh_db query results and saves the indexes associated to each id
                abbreviations_index = {vmh_db["results"][i]["abbreviation"]:i for i in range(len(vmh_db["results"]))}
                vmh_id = self.annotations["vmh"]
                met_vmh_index = abbreviations_index[vmh_id]
                #TODO: add exception handling
                inchiString = vmh_db["results"][met_vmh_index]["inchiString"]
                if inchiString is not None and inchiString != npNaN and inchiString != "":
                    self.inchi_string = inchiString
            elif "bigg.metabolite" in self.annotations.keys():
                met_db_id = self.annotations["bigg.metabolite"]
                bigg_db = pd.read_csv("../Databases/BiGG.tsv", sep='\t',index_col=0) 
                #TODO: continue BiGG handling and other databases...
            
                




        # series of if/else statements to handle all the types of annotations we have (can be handled by calling the functions of conversion between i.e. smiles to inchi string
            

        return None
    
    def compare(self, metabolite):

        return None
