# Porthmeus
# 08.03.24
from __future__ import annotations

from collections import defaultdict
from typing import Callable, Iterable, List, NewType, Sequence

import os
import sys
import warnings

import pandas as pd

from src.MeMoMetabolite import MeMoMetabolite
from src.annotation.annotateInchiRoutines import findOptimalInchi
from src.download_db import get_config, get_database_path

# HMDB0000972
AnnotationKey = NewType("AnnotationKey", str)
DBName = NewType("DBName", str)
DBKey = NewType("DBKey", str)

EntryAnnotationFunction = Callable[[AnnotationKey, pd.DataFrame, bool], tuple[dict, list, str]]


class AnnotationResult:
    def __init__(self, annotated_inchis: int, annotated_dbs: int, annotated_names: int):
        self.annotated_inchis: int = annotated_inchis
        self.annotated_dbs: int = annotated_dbs
        self.annotated_names: int = annotated_names
        self.annotated_total: int = annotated_inchis + annotated_dbs + annotated_names

    @classmethod
    def fromTuple(cls, annotationResult: tuple[int, int, int]) -> AnnotationResult:
        return cls(annotationResult[0], annotationResult[1], annotationResult[2])

    @classmethod
    def fromAnnotation(cls, annotationResult: AnnotationResult) -> AnnotationResult:
        return cls(
            annotationResult.annotated_inchis,
            annotationResult.annotated_dbs,
            annotationResult.annotated_names,
        )

    def __str__(self) -> str:
        return (
            f"Annotated inchis {self.annotated_inchis}, "
            f"annotated dbs {self.annotated_dbs}, annotated names {self.annotated_names}"
        )

    def __le__(self, other) -> bool:
        return (
            self.annotated_inchis <= other.annotated_inchis
            and self.annotated_dbs <= other.annotated_dbs
            and self.annotated_names <= other.annotated_names
        )

    def __gt__(self, other):
        return not self.__le__(other)

    def __add__(self, other):
        return AnnotationResult(
            self.annotated_inchis + other.annotated_inchis,
            self.annotated_dbs + other.annotated_dbs,
            self.annotated_names + other.annotated_names,
        )

    def __iter__(self):
        return iter([self.annotated_inchis, self.annotated_dbs, self.annotated_names])

    def __eq__(self, other):
        if not isinstance(other, AnnotationResult):
            return NotImplemented
        return (
            self.annotated_inchis == other.annotated_inchis
            and self.annotated_dbs == other.annotated_dbs
            and self.annotated_names == other.annotated_names
        )


####################################


def _resolve_database_filename(db_name: str) -> str:
    """Resolve a database key or filename to a filename under Databases/."""
    config = get_config()

    # db_name can be a database key (e.g. "BiGG") or already a filename
    if db_name in config["databases"]:
        db_cfg = config["databases"][db_name]
        raw_filename = db_cfg.get("file")
        reformat_filename = db_cfg.get("reformat")
        if not raw_filename and not reformat_filename:
            raise KeyError(f"No 'reformat' or 'file' entry configured for database '{db_name}'")

        # Prefer reformatted database if it exists on disk, otherwise fall back
        # to the raw download so annotation still works without reformatted DBs.
        if reformat_filename:
            reformat_path = os.path.join(get_database_path(), reformat_filename)
            if os.path.exists(reformat_path):
                return reformat_filename

        if raw_filename:
            return raw_filename

        # Last resort: return the reformat filename and let the caller handle
        # missing files according to allow_missing_dbs.
        return reformat_filename

    # Treat it as a direct filename
    return db_name


def load_database(
    db_name: str = "",
    allow_missing_dbs: bool = False,
    conversion_method: Callable[[str], pd.DataFrame] | None = None,
) -> pd.DataFrame:
    """
    Load the given database from the Databases folder.

    This function is backward-compatible with older call sites that passed a
    filename and a custom conversion method, and newer call sites that pass a
    configured database key (e.g. "BiGG").
    """
    try:
        filename = _resolve_database_filename(db_name)
        db_path = os.path.join(get_database_path(), filename)
        if conversion_method is not None:
            db = conversion_method(db_path)
        else:
            db = pd.read_csv(db_path, sep=",", header=0, dtype=str)
    except FileNotFoundError as e:
        if not allow_missing_dbs:
            raise e
        warnings.warn(str(e))
        db = pd.DataFrame()
    return db


def annotateEntry(entry: str, database: pd.DataFrame = pd.DataFrame()) -> tuple[dict, list, str]:
    """
    Retrieve annotations, names, and a best-effort InChI string for a database entry.

    The database is expected to follow the internal reformatted schema with columns
    like "id", "name", "DBs", and "inchi".
    """
    if database.empty:
        return dict(), list(), ""

    # Names
    if "name" in database.columns:
        names = database.loc[database["id"] == entry, "name"]
        all_names: list[str] = []
        for name in names:
            all_names.extend(name.split("_|_"))
    else:
        all_names = []

    # Cross-database annotations
    if "DBs" in database.columns:
        annos = database.loc[database["id"] == entry, "DBs"]
        all_annos: list[dict] = []
        for anno in annos:
            if isinstance(anno, str) and anno.startswith("{") and anno.endswith("}"):
                anno_dict = eval(anno)
            else:
                raise ValueError("Check your database files!")
            all_annos.append(anno_dict)
        merged_annos: dict[str, list] = defaultdict(list)
        for anno_dict in all_annos:
            for key, value in anno_dict.items():
                merged_annos[key].extend(value)
        merged_annos = dict(merged_annos)
    else:
        merged_annos = dict()

    # InChI
    opt_inchi = ""
    if "inchi" in database.columns:
        inchis = database.loc[database["id"] == entry, "inchi"].dropna()
        if len(inchis) != 0:
            opt_inchi = findOptimalInchi(inchis.tolist())

    return merged_annos, all_names, opt_inchi


def _unpack_annotation_function_result(result: Sequence) -> tuple[dict, list, str, str]:
    """Normalize annotation function return values to (annos, names, inchi, source)."""
    if len(result) == 4:
        annos, names, inchi, source = result
        return annos, names, inchi, source
    if len(result) == 3:
        annos, names, source = result
        return annos, names, "", source
    if len(result) == 2:
        annos, names = result
        return annos, names, "", ""
    raise ValueError(f"Unexpected annotation function result shape: {len(result)}")


def handleIDs(*args, **kwargs) -> AnnotationResult:
    """
    Backward-compatible ID-based annotation.

    Old style:
        handleIDs(db: pd.DataFrame, metabolites, db_key: str, annotation_function)

    New style:
        handleIDs(metabolites, db_name: str, allow_missing_dbs: bool = False)
    """
    # Old-style signature
    if args and isinstance(args[0], pd.DataFrame):
        db: pd.DataFrame = args[0]
        metabolites: List[MeMoMetabolite] = args[1]
        db_key: str = args[2]
        annotation_function: Callable = args[3]

        new_annos = 0
        new_names = 0
        new_inchis = 0

        for met in metabolites:
            if met._id is None or db.empty or db_key not in db.columns:
                continue
            if not any(db[db_key] == met._id):
                continue

            annos, names, inchi, source = _unpack_annotation_function_result(
                annotation_function(met._id, db)
            )
            source = source or "id"

            if names:
                new_names += met.add_names(names, source=source)
            if annos:
                new_annos += met.add_annotations(annos, source=source)
            if inchi:
                new_inchis += met.add_inchi_string(inchi, source)

        return AnnotationResult(new_inchis, new_annos, new_names)

    # New-style signature
    metabolites: List[MeMoMetabolite] = args[0]
    db_name: str = args[1]
    allow_missing_dbs: bool = kwargs.get("allow_missing_dbs", False)

    db = load_database(db_name, allow_missing_dbs)
    if db.empty or "id" not in db.columns:
        return AnnotationResult(0, 0, 0)

    new_annos = 0
    new_names = 0
    new_inchis = 0
    source = db_name

    for met in metabolites:
        if met._id is None:
            continue
        if not any(db["id"] == met._id):
            continue

        annos, names, inchi = annotateEntry(met._id, db)

        if names:
            new_names += met.add_names(names, source)
        if annos:
            new_annos += met.add_annotations(annos, source)
        if inchi:
            new_inchis += met.add_inchi_string(inchi, source)

    return AnnotationResult(new_inchis, new_annos, new_names)


def handleMetabolites(*args, **kwargs) -> AnnotationResult:
    """
    Backward-compatible annotation via metabolite annotation slots.

    Old style:
        handleMetabolites(db: pd.DataFrame, metabolites, db_key: str, annotation_function)

    New style:
        handleMetabolites(metabolites, db_name: str, allow_missing_dbs: bool = False)
    """
    # Old-style signature
    if args and isinstance(args[0], pd.DataFrame):
        db: pd.DataFrame = args[0]
        metabolites: List[MeMoMetabolite] = args[1]
        db_key: str = args[2]
        annotation_function: Callable = args[3]

        new_annos_added = 0
        new_names_added = 0
        new_inchis_added = 0

        for met in metabolites:
            if db.empty or db_key not in met.annotations:
                continue

            new_met_anno: dict[str, list] = {}
            new_names: list[str] = []
            new_inchi = ""
            source = "annotation"

            for entry in met.annotations[db_key]:
                annos, names, inchi, src = _unpack_annotation_function_result(
                    annotation_function(entry, db)
                )
                source = src or source

                for key, value in annos.items():
                    if key in new_met_anno:
                        new_met_anno[key].extend(value)
                    else:
                        new_met_anno[key] = value
                new_names.extend(names)
                if inchi:
                    new_inchi = findOptimalInchi([new_inchi, inchi])

            if new_names:
                new_names_added += met.add_names(new_names, source)
            if new_met_anno:
                new_annos_added += met.add_annotations(new_met_anno, source)
            if new_inchi:
                new_inchis_added += met.add_inchi_string(new_inchi, source)

        return AnnotationResult(new_inchis_added, new_annos_added, new_names_added)

    # New-style signature
    metabolites: List[MeMoMetabolite] = args[0]
    db_name: str = args[1]
    allow_missing_dbs: bool = kwargs.get("allow_missing_dbs", False)

    # Conversion from internal database name to identifiers.org namespace
    db_keys = {
        "BiGG": "bigg.metabolite",
        "VMH": "vmhmetabolite",
        "HMDB": "HMDB",
        "ModelSeed": "seed.compound",
        "ChEBI": "chebi",
        "gapseq": "seed.compound",
    }
    db_key = db_keys.get(db_name)
    if db_key is None:
        return AnnotationResult(0, 0, 0)

    db = load_database(db_name, allow_missing_dbs)
    if db.empty:
        return AnnotationResult(0, 0, 0)

    new_annos_added = 0
    new_names_added = 0
    new_inchis_added = 0
    source = db_name

    for met in metabolites:
        if db_key not in met.annotations:
            continue

        new_met_anno: dict[str, list] = {}
        new_names: list[str] = []
        new_inchi = ""

        for entry in met.annotations[db_key]:
            annos, names, inchi = annotateEntry(entry, db)
            for key, value in annos.items():
                if key in new_met_anno:
                    new_met_anno[key].extend(value)
                else:
                    new_met_anno[key] = value
            new_names.extend(names)
            if inchi:
                new_inchi = findOptimalInchi([new_inchi, inchi])

        if new_names:
            new_names_added += met.add_names(new_names, source)
        if new_met_anno:
            new_annos_added += met.add_annotations(new_met_anno, source)
        if new_inchi:
            new_inchis_added += met.add_inchi_string(new_inchi, source)

    return AnnotationResult(new_inchis_added, new_annos_added, new_names_added)
