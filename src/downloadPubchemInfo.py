# Porthmeus
# 15.02.23


"""
Get a list of metabolite names, search pubchem for it and retrieve an InChI/Smile string from it.
"""


import pubchempy as pcp
import pandas as pd
import numpy as np
import time
import sys
import warnings
import urllib.error

# global variable for the atoms list
atoms_list = [
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
    "R",
]


def downloadElements() -> pd.DataFrame:
    """
    Download the elements table from pubchem and return it as a data frame
    """

    # download the element list from pubchem
    try:
        elements = pd.read_csv(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/periodictable/CSV?response_type=save&response_basename=PubChemElements_all"
        )
    except FileNotFoundError as e:
        warnings.warn(e)
        elements = pd.DataFrame()
    return elements


def looksLikeSumFormula(string: str, atoms: list[str] = atoms_list) -> bool:
    """
    Takes a string and tries to determine if it is likely a sum formula
    @atoms: a list of all atom symbols from the elements table
    """

    # add the R for undifined residues
    atoms.append("R")

    # take the string an split it by numbers and capital letters and check whether the substrings are atoms
    isAtom = []
    atom = ""
    for letter in string:
        if not letter.isdigit():
            if letter.isupper():
                if atom != "":
                    if len(atom) < 3:
                        isAtom.append(int(atom in atoms))
                    else:
                        isAtom.append(0)
                atom = letter
            else:
                atom = atom + letter
    # check for the last entry in the list
    if atom != "":
        if len(atom) < 3:
            isAtom.append(int(atom in atoms))
        else:
            isAtom.append(0)

    # check if at least 80% of the strings found are atoms than report that it is likely a sum formula, otherwise report false
    if sum(isAtom) / len(isAtom) >= 0.8:
        return True
    else:
        return False


def removeSumFormulaFromName(name: str, atoms: list[str] = atoms_list) -> str:
    """
    Takes a name of a metabolite and removes sum formulas, if they should be given for some reasons
    """
    name_l = [x for x in name.split() if not looksLikeSumFormula(x, atoms)]
    return " ".join(name_l)


def downloadCompInfo(metabolite: str) -> pd.DataFrame:
    try:
        comp = pcp.get_compounds(metabolite, "name", as_dataframe=True)
    except KeyError:
        # add here more manipulations to the metabolite name if there is a problem with it
        metabolite_rm = removeSumFormulaFromName(metabolite)
        if metabolite_rm == "":
            metabolite_rm = metabolite.split()[0]
        sys.stdout.write(
            "\r Trying to remove sum formulas from name. "
            + metabolite
            + " --> "
            + metabolite_rm
        )
        sys.stdout.flush()
        try:
            comp = pcp.get_compounds(metabolite_rm, "name", as_dataframe=True)
        except KeyError:
            try:
                # if that does not help, see if a substance exists in the database
                comp = pcp.get_substances(metabolite_rm, "name", as_dataframe=True)
            except:
                # if the problem persists return an empty data frame
                comp = pd.DataFrame()
    return comp


def getPubchemInfo(metabolites: list[str]) -> pd.DataFrame:
    """
    Go through a list of metabolites and return the info from a pubchem search of them as a pandas dataframe
    """
    # for removing sum formulas from names load the atoms list from pubchem elements table
    comp_all = pd.DataFrame()
    n = len(metabolites)
    try:
        for i, met_n in enumerate(metabolites):
            comp = downloadCompInfo(met_n)
            comp = comp.assign(Name=met_n)
            comp_all = pd.concat([comp_all, comp])
            time.sleep(1)
            sys.stdout.write(
                "\r" + str(round(((i + 1) / n) * 100, ndigits=2)) + "% done" + " " * 80
            )
            sys.stdout.flush()
    except urllib.error.URLError as e:
        print("Could not establish connection to PubMed")
        return pd.DataFrame()

    return comp_all
