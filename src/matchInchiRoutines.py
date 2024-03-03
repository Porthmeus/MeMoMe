# Porthmeus
# 15.03.23

# Some functions to compare InChi keys/strings and tell if it is the same chemical structure

import rdkit
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.rdchem import Mol


def compareInchiByFingerprint(met1: str, met2: str, method: str = "Dice", verbose = False) -> float:
    ''' Simple function to compare based on fingerprints
    @met1 - str of the first inchi string
    @met2 - str of the second inchi string
    @method - str defining the method to use. One of 'Tanimoto', 'Dice', 'Cosine', 'Sokal', 'Russel', 'Kulczynski', 'McConnaughey'.'''
    
    # turn off chattiness of rdkit
    if verbose == False:
        RDLogger.DisableLog("rdApp.*")

    # correct and check method string if implemented in RDKit
    # method = method[0].upper() + method[1:].lower()
    impl_methods = ['Tanimoto', 'Dice', 'Cosine', 'Sokal', 'Russel', 'Kulczynski',
                    'McConnaughey']  # "Tversky" is noted in the documentation, but throws an error
    if method not in impl_methods:
        raise ValueError(method + " not implemented. Choose one of: " + ", ".join(impl_methods))

    # do the conversion
    m1 = Chem.MolFromInchi(met1)
    m2 = Chem.MolFromInchi(met2)
    fp1 = Chem.RDKFingerprint(m1)
    fp2 = Chem.RDKFingerprint(m2)
    sim = DataStructs.FingerprintSimilarity(fp1, fp2, metric=eval("DataStructs." + method + "Similarity"))
    return (sim)

def compareInchiByFingerprint0(fp1, fp2, method: str = "Dice", verbose:bool = False) -> float:
    ''' Simple function to compare based on fingerprints
    @met1 - str of the first inchi string
    @met2 - str of the second inchi string
    @method - str defining the method to use. One of 'Tanimoto', 'Dice', 'Cosine', 'Sokal', 'Russel', 'Kulczynski', 'McConnaughey'.'''
    
    # turn off chattiness of rdkit
    if verbose == False:
        RDLogger.DisableLog("rdApp.*")

    # correct and check method string if implemented in RDKit
    # method = method[0].upper() + method[1:].lower()
    impl_methods = ['Tanimoto', 'Dice', 'Cosine', 'Sokal', 'Russel', 'Kulczynski',
                    'McConnaughey']  # "Tversky" is noted in the documentation, but throws an error
    if method not in impl_methods:
        raise ValueError(method + " not implemented. Choose one of: " + ", ".join(impl_methods))

    # do the conversion
    sim = DataStructs.FingerprintSimilarity(fp1, fp2, metric=eval("DataStructs." + method + "Similarity"))
    return (sim)

def compareInchiByCSmile(met1: str, met2: str, verbose = False) -> bool:
    '''Converts Inchi to canonical smile and compares the strings'''

    # turn off chattiness of rdkit
    if verbose == False:
        RDLogger.DisableLog("rdApp.*")
    # do the conversion
    m1 = Chem.MolFromInchi(met1)
    m2 = Chem.MolFromInchi(met2)
    sm1 = Chem.MolToSmiles(m1)
    sm2 = Chem.MolToSmiles(m2)
    return (sm1 == sm2)


def compareInchiByKey(met1: str, met2: str, verbose = False) -> bool:
    '''Converts Inchi to Inchi-key and compares the hashes'''

    # turn off chattiness of rdkit
    if verbose == False:
        RDLogger.DisableLog("rdApp.*")

    # do the conversion
    m1 = Chem.MolFromInchi(met1)
    m2 = Chem.MolFromInchi(met2)
    key1 = Chem.MolToInchiKey(m1)
    key2 = Chem.MolToInchiKey(m2)
    return (key1 == key2)


def compareInchiByString(met1: str, met2: str, verbose = False) -> bool:
    '''Converts Inchi to canonical Inchi and compares the strings'''

    # turn off chattiness of rdkit
    if verbose == False:
        RDLogger.DisableLog("rdApp.*")

    # do the conversion
    m1 = Chem.MolFromInchi(met1)
    m2 = Chem.MolFromInchi(met2)
    key1 = Chem.MolToInchi(m1)
    key2 = Chem.MolToInchi(m2)
    return (key1 == key2)

def compareInchiByString0(m1:rdkit.Chem.rdchem.Mol, m2:rdkit.Chem.rdchem.Mol, verbose:bool = False) -> bool:
    '''Converts Inchi to canonical Inchi and compares the strings'''

    # turn off chattiness of rdkit
    if verbose == False:
        RDLogger.DisableLog("rdApp.*")

    # do the conversion
    key1 = Chem.MolToInchi(m1)
    key2 = Chem.MolToInchi(m2)
    return (key1 == key2)

def compareInchiByStringNoHydrogen(met1: str, met2: str, verbose = False) -> bool:
    '''Converts Inchi to canonical InChI after removal of Hydrogens'''

    # turn off chattiness of rdkit
    if verbose == False:
        RDLogger.DisableLog("rdApp.*")

    # do conversion
    m1 = Chem.MolFromInchi(met1)
    m2 = Chem.MolFromInchi(met2)
    # remove Hs
    m1 = Chem.RemoveAllHs(m1)
    m2 = Chem.RemoveAllHs(m2)
    # compare string
    key1 = Chem.MolToInchi(m1)
    key2 = Chem.MolToInchi(m2)
    return key1 == key2


def compareInchiByStereoIsomer(met1: str, met2: str, verbose = False) -> bool:
    '''Enumerates all stereoisomers of both metabolites and checks if there is a common set between both metabolites'''

    # turn off chattiness of rdkit
    if verbose == False:
        RDLogger.DisableLog("rdApp.*")

    # do conversion
    m1 = Chem.MolFromInchi(met1)
    m2 = Chem.MolFromInchi(met2)
    # find also stereoisomers in defined InChI
    opts = StereoEnumerationOptions(onlyUnassigned=False)
    stereo1 = [Chem.MolToInchi(x) for x in EnumerateStereoisomers(m1, options=opts)]
    stereo2 = [Chem.MolToInchi(x) for x in EnumerateStereoisomers(m2, options=opts)]
    intersect = list(set(stereo1).intersection(stereo2))
    return len(intersect) > 0

def compareInchiByStereoIsomer0(m1:rdkit.Chem.rdchem.Mol, m2:rdkit.Chem.rdchem.Mol, verbose:bool = False) -> bool:
    '''Enumerates all stereoisomers of both metabolites and checks if there is a common set between both metabolites'''

    # turn off chattiness of rdkit
    if verbose == False:
        RDLogger.DisableLog("rdApp.*")
    
    # sanitize molecules
    Chem.SanitizeMol(m1)
    Chem.SanitizeMol(m2)
     
    # find also stereoisomers in defined InChI
    opts = StereoEnumerationOptions(onlyUnassigned=False)
    stereo1 = [Chem.MolToInchi(x) for x in EnumerateStereoisomers(m1, options=opts)]
    stereo2 = [Chem.MolToInchi(x) for x in EnumerateStereoisomers(m2, options=opts)]
    intersect = list(set(stereo1).intersection(stereo2))
    return len(intersect) > 0


def NeutraliseCharges(mol: rdkit.Chem.rdchem.Mol, verbose = False) -> rdkit.Chem.rdchem.Mol:
    """ Takes a molecule and returns the neutralized version of it + an indicator, what the charge difference was."""
    
    # turn off chattiness of rdkit
    if verbose == False:
        RDLogger.DisableLog("rdApp.*")

    # """ adapted from Hans de Winter - https://rdkit.readthedocs.io/en/latest/Cookbook.html """
    # initialize the patterns
    patts = (
        # Imidazoles
        ('[n+;H]', 'n'),
        # Amines
        ('[N+;!H0]', 'N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]', 'O'),
        # Thiols
        ('[S-;X1]', 'S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]', 'N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]', 'N'),
        # Tetrazoles
        ('[n-]', '[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]', 'S'),
        # Amides
        ('[$([N-]C=O)]', 'N'),
    )
    reactions = [(Chem.MolFromSmarts(x), Chem.MolFromSmiles(y, False)) for x, y in patts]

    for i, (reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    return (mol)
