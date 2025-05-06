"""
This file contains function to annotated MeMoMetabolites in a bulk manner. This means a list of MeMoMetabolites is parsed and the annotation takes place on the metabolites in that list. This is mainly to use within the MeMoModel.annotate() function
"""

import os
import json
import pandas as pd
import warnings
from io import StringIO
from src.download_db import get_config, get_database_path
from src.MeMoMetabolite import MeMoMetabolite
from src.annotation.annotateAux import (
    AnnotationResult,
    load_database,
    handleIDs,
    handleMetabolites,
    DBName,
    DBKey,
    AnnotationKey,
)
from src.annotation.annotateALL import *
from typing import Optional, Tuple


def handle_vmh_entries(vmh: pd.DataFrame, entry: str) -> Tuple[dict, list[str], str]:
    keys = [
        "biggId",
        "keggId",
        "cheBlId",
        "inchiString",
        "inchiKey",
        "smile",
        "hmdb",
        "metanetx",
        "seed",
        "biocyc",
    ]
    vmh = vmh.loc[vmh["abbreviation"] == entry,]
    fullName = list(set(vmh["fullName"]))
    vmh = vmh[keys]
    annotations = dict()
    if len(vmh) > 1:
        raise Exception("More than one value for key: " + entry)

    for index, row in vmh.iterrows():
        for key, value in row.items():
            if (not pd.isna(value)) and (len(str(value)) > 0):
                # If key does not exist, create it and append the value
                annotations.setdefault(key, []).append(value)

    return (annotations, fullName, "VMH")


def __json_to_tsv(path: str) -> pd.DataFrame:
    with open(path, "r") as f:
        vmh = json.load(f)["results"]
        return pd.DataFrame(vmh)


def annotateVMH_entry(
    entry: AnnotationKey,
    database: pd.DataFrame = pd.DataFrame(),
    allow_missing_dbs: bool = False,
) -> tuple[dict, list, str]:
    """
    A small helper function to avoid redundant code.
    Uses a VMH identifier and annotates it with the identifiers.org entries.
    database - is either empty (default) or a dictionary with the data from the VMH homepage for the metabolites. If it is empty, the function will try to load the database from the config file.
    Returns a tuple containing a dictionary for the extracted annotations and a list of new names for the metabolite.
    """
    return annotateDB_entry(
        entry,
        DBName("VMH"),
        __json_to_tsv,
        handle_vmh_entries,
        database,
        allow_missing_dbs,
    )


def annotateVMH(
    metabolites: list[MeMoMetabolite], allow_missing_dbs: bool = False
) -> AnnotationResult:
    """
    Annotate a list of metabolites with the entries from VMH. Look for VMH ids in the annotation slot of the metabolites and if one is found use these.
    The function will directly add the annotation and names to the MeMoMetabolite object.
    Return value: a tuple of 3 denoting the number of changes metabolites for the following values: [inchi_strings, annotations, names]
    """

    vmh = load_database(
        get_config()["databases"]["VMH"]["file"], allow_missing_dbs, __json_to_tsv
    )

    if vmh.empty:
        return AnnotationResult(0, 0, 0)

    return handleMetabolites(vmh, metabolites, "vmhmetabolite", annotateVMH_entry)


def annotateVMH_id(
    metabolites: list[MeMoMetabolite], allow_missing_dbs: bool = False
) -> AnnotationResult:
    """
    Annotate a list of metabolites with the entries from VMH. Look for VMH ids in the metabolite._id slot and if one is found use these.
    """

    # check if the database was given, if not, try to load it
    vmh = load_database(
        get_config()["databases"]["VMH"]["file"], allow_missing_dbs, __json_to_tsv
    )

    if vmh.empty:
        return AnnotationResult(0, 0, 0)

    return handleIDs(vmh, metabolites, "abbreviation", annotateVMH_entry)
