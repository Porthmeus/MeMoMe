Class:Metabolite
values:
    - ID (original from model)
    - Model (which model does it originate from)
    - Names [...]
    - InChI-string
    - Charge
    - Databases {DB:[ids]}
    - Attributes {...}
methods:
    - annotate(self): try to get the InChI with different methods
    - compare(self, Metabolite): -> [bool, score, method] compare itself to other metabolite classes
         - compare by InChI
         - compare by DB ids
         - compare by fuzzy name matching
