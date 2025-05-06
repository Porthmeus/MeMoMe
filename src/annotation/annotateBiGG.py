# Porthmeus
# 16.02.24

"""
This file contains function to annotated MeMoMetabolites in a bulk manner. This means a list of MeMoMetabolites is parsed and the annotation takes place on the metabolites in that list. This is mainly to use within the MeMoModel.annotate() function
"""

import os
import warnings
import pandas as pd
from src.MeMoMetabolite import MeMoMetabolite
from src.download_db import get_config, get_database_path
from src.parseMetaboliteInfos import getAnnoFromIdentifierURL
from src.annotation.annotateAux import (
    AnnotationResult,
    load_database,
    handleIDs,
    handleMetabolites,
    AnnotationKey,
)
from typing import Dict, List


def handle_bigg_entries(urls: pd.Series) -> Dict[str, List[str]]:
    """
    Internal method to handle the annotation of bigg entrires, that maps biggId to identifier.org urls.
    This method handles the db/id extraction from those urls
    urls: pd.Series: This should be a series with just one column and multiple rows.
      Each cell should be a list of urls
    """
    annotations = dict()
    if len(urls) > 0:
        # Urls is now a df that has one column and multiple rows
        # Each row is a bag of urls
        # To convert them to a list we join them by semicolon and again split them
        # that results is a list that contains all the urls
        urls = ";".join(list(urls))
        urls = urls.split(";")
        for url in urls:
            # Src is a db, e.g. hmdb or seed
            # val is a identifier in the corresponding db
            src, val = getAnnoFromIdentifierURL(url)
            if src in annotations.keys():
                annotations[src].append(val)
            elif src != None:
                annotations[src] = [val]

    return annotations


def annotateBiGG_entry(
    entry: AnnotationKey,
    database: pd.DataFrame = pd.DataFrame(),
    allow_missing_dbs: bool = False,
) -> tuple[dict, list, str]:
    """
    A small helper function to avoid redundant code
    Uses a BiGG identifiers and annotates it with the identifiers.org entries.
    database - is either empty (default) or a pandas data.frame with the data from the BiGG model homepage for the metabolites. If it is empty, the function will try to load the database from the config file
    Returns a tuple containing a dictionary for the extracted annotations and a list of new names for the metabolite
    """
    # check if the database was given, if not, try to load it
    if len(database) == 0:
        bigg = load_database(
            get_config()["databases"]["BiGG"]["file"],
            allow_missing_dbs,
            lambda path: pd.read_csv(path, sep="\t"),
        )
    else:
        bigg = database

    if bigg.empty:
        return dict(), list(), "bigg"

    urls = bigg.loc[bigg["universal_bigg_id"] == entry, "database_links"]
    urls = urls.loc[~pd.isna(urls)]

    names = list(pd.unique(bigg.loc[bigg["universal_bigg_id"] == entry, "name"]))
    annotations = handle_bigg_entries(urls)
    return annotations, names, "bigg"


def annotateBiGG(
    metabolites: list[MeMoMetabolite], allow_missing_dbs: bool = False
) -> AnnotationResult:
    """
    Annotate a list of metabolites with the entries from BiGG. Look for BiGG ids in the annotation slot of the metabolites and if one is found use these. Since BiGG does not provide any InChI strings, this will not results in any new annotated Inchi strings, but it will increase the number of entries in the metabolite annotation.
    The function will directly add the annotation and names to the MeMoMetabolite object.
    Return value: a tuple of 3 denoting the number of changes metabolites for the following values: [inchi_strings, annotations, names]
    """

    # load the database
    bigg = load_database(
        get_config()["databases"]["BiGG"]["file"],
        allow_missing_dbs,
        lambda path: pd.read_csv(path, sep="\t"),
    )
    if bigg.empty:
        return AnnotationResult(0, 0, 0)

    return handleMetabolites(bigg, metabolites, "bigg.metabolite", annotateBiGG_entry)


def annotateBiGG_id(
    metabolites: list[MeMoMetabolite], allow_missing_dbs: bool = False
) -> AnnotationResult:
    """
    Annotate a list of metabolites with the entries from BiGG. Look for BiGG ids in the metabolite._id slot and if one is found use these. Since BiGG does not provide any InChI strings, this will not results in any new annotated Inchi strings, but it will increase the number of entries in the metabolite annotation.
    """

    bigg = load_database(
        get_config()["databases"]["BiGG"]["file"],
        allow_missing_dbs,
        lambda path: pd.read_csv(path, sep="\t"),
    )
    if bigg.empty:
        return AnnotationResult(0, 0, 0)

    return handleIDs(bigg, metabolites, "universal_bigg_id", annotateBiGG_entry)
