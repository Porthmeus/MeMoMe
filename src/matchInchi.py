# Porthmeus
# 15.03.23

# Some functions to compare InChi keys/strings and tell if it is the same chemical structure

from rdkit import Chem,DataStructs
from src.downloadPubchemInfo import getPubchemInfo



def compareInchiByFingerprint(met1:str, met2:str, method:str = "Dice") -> float:
    ''' Simple function to compare based on fingerprints
    @met1 - str of the first inchi string
    @met2 - str of the second inchi string
    @method - str defining the method to use. One of 'Tanimoto', 'Dice', 'Cosine', 'Sokal', 'Russel', 'Kulczynski', 'McConnaughey'.'''
    
    # correct and check method string if implemented in RDKit
    method = method[0].upper() + method[1:].lower()
    impl_methods = ['Tanimoto', 'Dice', 'Cosine', 'Sokal', 'Russel', 'Kulczynski', 'McConnaughey'] # "Tversky" is noted in the documentation, but throws an error
    if method not in impl_methods:
        raise ValueError(method+ " not implemented. Choose one of: " + ", ".join(impl_methods))

    # do the conversion
    m1 = Chem.MolFromInchi(met1)
    m2 = Chem.MolFromInchi(met2)
    fp1 = Chem.RDKFingerprint(m1)
    fp2 = Chem.RDKFingerprint(m2)
    sim = DataStructs.FingerprintSimilarity(fp1,fp2, metric = eval("DataStructs."+method+"Similarity"))
    return(sim)

def compareInchiByCSmile(met1:str, met2:str) -> bool:
    '''Converts Inchi to canonical smile and compares the strings'''
    # do the conversion
    m1 = Chem.MolFromInchi(met1)
    m2 = Chem.MolFromInchi(met2)
    sm1 = Chem.MolToSmiles(m1)
    sm2 = Chem.MolToSmiles(m2)
    return(sm1 == sm2)


def compareInchiByKey(met1:str, met2:str) -> bool:
    '''Converts Inchi to Inchi-key and compares the hashes'''
    # do the conversion
    m1 = Chem.MolFromInchi(met1)
    m2 = Chem.MolFromInchi(met2)
    key1 = Chem.MolToInchiKey(m1)
    key2 = Chem.MolToInchiKey(m2)
    return(key1 == key2)

def compareInchiByString(met1:str, met2:str) -> bool:
    '''Converts Inchi to canonical Inchi and compares the strings'''
    # do the conversion
    m1 = Chem.MolFromInchi(met1)
    m2 = Chem.MolFromInchi(met2)
    key1 = Chem.MolToInchi(m1)
    key2 = Chem.MolToInchi(m2)
    return(key1 == key2)
