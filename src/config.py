from enum import Enum
from typing import Dict


class Database(Enum):
    ModelSeed = "modelSeed.tsv"
    BiGG = "BiGG.tsv"
    VMH = "vmh.json"


database_links: Dict[Database, str] = {
    Database.ModelSeed: "https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/compounds.tsv",
    Database.BiGG: "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt",
    Database.VMH: "https://www.vmh.life/_api/metabolites/?format=jsonp&page_size=9999999&_dc=1670601370202&page=1&start=0&limit=9999999&callback=Ext.data.JsonP.callback19"
}
