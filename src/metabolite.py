#Class:Metabolite
#values:
#    - Attributes {...}
#methods:
#    - annotate(self): try to get the InChI with different methods
#    - compare(self, Metabolite): -> [bool, score, method] compare itself to other metabolite classes
#         - compare by InChI
#         - compare by DB ids
#         - compare by fuzzy name matching

            
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
    databases: dict()
       A dictionary whose keys are the names of the databases, and the items are the database-specific metabolite ids
    """

    def __init__(
        self,
        orig_id:str = None,
        model:str = None,
        names:list = [],
        inchi_string:str = None,
        formula:str = None,
        charge:float = None,
        databases:dict = None,
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
        databases: dict()
           A dictionary whose keys are the names of the databases, and the items are the database-specific metabolite ids
        """
        self.orig_id = orig_id
        self.model = model
        self.names = names
        self.formula = formula
        self.inchi_string = inchi_string
        self.charge = charge
        self.databases = databases

        
    def annotate(self):
        if "bigg.metabolite" in self.databases.keys():
            met_db_id = self.databases["bigg.metabolite"]
            bigg_db = pd.read_csv("../Databases/BiGG.tsv", sep='\t',index_col=0)
            

        return None
    
    def compare(self, metabolite):
        return None
