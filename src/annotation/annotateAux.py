# Porthmeus
# 08.03.24
from __future__ import annotations
from collections import defaultdict
import pandas as pd
import warnings
import os
import sys
from typing import Callable, List, NewType
from src.download_db import get_config, get_database_path
from src.MeMoMetabolite import MeMoMetabolite
from src.annotation.annotateInchiRoutines import findOptimalInchi

# HMDB0000972
AnnotationKey = NewType('AnnotationKey', str)
# Database name, corresponds to one of the keys in config.yaml
DBName = NewType('DBName', str)
# The dbs have have certain columns that correspond to their main identifier, e.g. vmh has vmhmetabolite, HMDB has accession and so one. This is what we call a database key. Aka the MAIN column which we use for entry lookup
DBKey = NewType('DBKey', str)

EntryAnnotationFunction = Callable[[AnnotationKey, pd.DataFrame, bool], tuple[dict, list, str]]

def annotateEntry(entry,
                  database:pd.DataFrame = pd.DataFrame()) -> tuple[dict, list, str, str | None]:
    """Takes an entry (id) for the database (db_name) and retrieves the annotations from the database.
    Return: A tuple with a dictionary containing the crosslinks to other databases, a list with alternative trivial names for the metabolite, the InChI string and the string from which database the annotation was obtained"""
    if database.empty:
        return dict(), list(), "", None
    
    # get the names of the entry
    if "name" in database.columns:
        names = database.loc[database["id"] == entry, "name"]
        all_names = list()
        for name in names:
            if not pd.isna(name):
                all_names.extend(name.split("_|_"))
    else:
        all_names = []
    all_names = sorted(set(all_names))
    
    # get annotations
    if "DBs" in database.columns:
        annos = database.loc[database["id"] == entry, "DBs"]
        all_annos = []
        for anno in annos:
            if anno.startswith("{") and anno.endswith("}"):
                anno = eval(anno)
            else:
                raise ValueError("Check your database files!")
            all_annos.append(anno)
        merged_annos = defaultdict(list)
        for d in all_annos:
            for k,v in d.items():
                merged_annos[k].extend(v)
        merged_annos = dict(merged_annos)
    else:
        merged_annos = dict()
    # get the inchi string
    if "inchi" in database.columns:
        inchis = database.loc[database["id"] == entry, "inchi"]
        inchis = inchis.dropna()
        if len(inchis) != 0:
            opt_inchi = findOptimalInchi(inchis.tolist())
        else:
            opt_inchi = ""
    
    if "formula" in database.columns:
        formula_series = database.loc[database["id"] == entry, "formula"]
        formula_series = formula_series.dropna()

        if not formula_series.empty:
          formula = str(formula_series.iloc[0])
        else:
          formula = None# Fallback if it was empty or NaN

    else:
      formula = None

    return merged_annos, all_names, opt_inchi, formula

class AnnotationResult():
  def __init__(self, annotated_inchis: int, annotated_dbs: int, annotated_names: int, annotated_formulas: int):
    self.annotated_inchis: int =  annotated_inchis
    self.annotated_dbs: int = annotated_dbs
    self.annotated_names: int = annotated_names
    self.annotated_formulas: int = annotated_formulas 
    self.annotated_total: int = annotated_inchis + annotated_dbs + annotated_names + annotated_formulas
  
  @classmethod
  def fromTuple(cls, annotationResult: tuple[int, int, int, int]) -> AnnotationResult:
    return cls(annotationResult[0], annotationResult[1], annotationResult[2], annotationResult[3])
  
  @classmethod
  def fromAnnotation(cls, annotationResult: AnnotationResult) -> AnnotationResult:
    return cls(annotationResult.annotated_inchis, annotationResult.annotated_dbs, annotationResult.annotated_names, annotationResult.annotated_formulas)

  def __str__(self) -> str:
    return f"Annotated inchis {self.annotated_inchis}, annotated dbs {self.annotated_dbs}, annotated names {self.annotated_names} Annotated formulas {self.annotated_formulas}"

  def __repr__(self) -> str:
    return f"Annotated inchis {self.annotated_inchis}, annotated dbs {self.annotated_dbs}, annotated names {self.annotated_names} Annotated formulas {self.annotated_formulas}"

  def __le__(self, other) -> bool:
    return self.annotated_inchis <= other.annotated_inchis and self.annotated_dbs <= other.annotated_dbs and self.annotated_names <= other.annotated_names and self.annotated_formulas <= other.annotated_formulas

  def __gt__(self, other):
    return not (self.__le__(other))

  def __add__(self, other):
    return AnnotationResult(self.annotated_inchis + other.annotated_inchis, self.annotated_dbs + other.annotated_dbs, self.annotated_names + other.annotated_names, self.annotated_formulas + other.annotated_formulas)

  def __iter__(self):
      return iter([self.annotated_inchis, self.annotated_dbs, self.annotated_names, self.annotated_formulas])

  def __eq__(self, other): 
      if not isinstance(other, AnnotationResult):
          # don't attempt to compare against unrelated types
          return NotImplemented

      return self.annotated_inchis == other.annotated_inchis and self.annotated_dbs == other.annotated_dbs and self.annotated_names == other.annotated_names, self.annotated_formulas == other.annotated_formulas
####################################


def load_database(db_name: str = "", allow_missing_dbs: bool = False) -> pd.DataFrame:
  """
  Load the given database. The file should in the projects root /Database folder.
  """
  # load the database
  try:
    db_path =  os.path.join(get_database_path(), get_config()["databases"][db_name]["reformat"])
    db = pd.read_csv(db_path, sep = ",", header = 0, dtype = str)
  except FileNotFoundError as e:
    # Rethrow exception if we don't allow missing dbs
    if allow_missing_dbs == False:
      raise e
    warnings.warn(str(e))
    db = pd.DataFrame()
  return(db)

def handleIDs(metabolites: List[MeMoMetabolite], db_name:DBName,  allow_missing_dbs: bool = False) -> AnnotationResult:
  """
  Checks for each metabolite if the metabolite id can be found in the column `db_key` can be found in `db`. 
  db: a dataframe with columns `db_key`. This column will be comparted to the metabolite id (met._id)
  metabolites: A list of metabolites that will be checked
  db_key: The column in the db dataframe
  annotation_function: Defines how to get an entry from the db which the given met._id (Check annotateVMH/BiGG for example usages.
  """
  db = load_database(db_name, allow_missing_dbs)
  new_annos = 0
  new_names = 0
  new_inchis = 0
  new_formulas = 0
  source = db_name
  for met in metabolites:
    if any(db["id"]==met._id):
      if met._id is None: raise Exception("met._id is None")
      new_met_anno_entry, new_names_entry, new_inchi_entry, new_formula = annotateEntry(met._id,  db)
      # add names
      if len(new_names_entry) >0:
          x = met.add_names(new_names_entry, source)
          new_names = new_names + x

      # add the annotations to the slot in the metabolites
      if len(new_met_anno_entry) > 0:
           x = met.add_annotations(new_met_anno_entry, source)
           new_annos = new_annos + x 

      if new_inchi_entry != "":
           x = met.add_inchi_string(new_inchi_entry, source)
           new_inchis = new_inchis + x
      
      if new_formula is not None: 
           x = met.set_formula(new_formula, source)
           new_formulas = new_formulas + x
  # return the number of metabolites which got newly annotated with inchis,
  # annotations and names
  anno_result = AnnotationResult(new_inchis, new_annos, new_names, new_formulas)
  return anno_result


def handleMetabolites(metabolites: List[MeMoMetabolite],  db_name:DBName, allow_missing_dbs: bool = False) -> AnnotationResult:
    """Checks the annotation dictionary entries and use them for further annotation of the metabolites"""
    # conversion from internal database annotation to identifier.org annotation
    # TODO save this as yaml to enable user annotation tables
    db_keys = {"BiGG" : "bigg.metabolite",
               "VMH" : "vmhmetabolite",
               "HMDB" : "HMDB",
               "ModelSeed" : "seed.compound",
               "ChEBI" : "chebi"}
    db_key = db_keys[db_name]
    db = load_database(db_name, allow_missing_dbs)
    new_annos_added = 0
    new_names_added = 0
    new_inchis_added = 0
    new_formulas_added = 0
    # go through the metabolites and check if there is data whigh can be added
    source = db_name
    for met in metabolites:
        new_met_anno = dict()
        new_names = list() 
        new_inchi = ""
        new_formula = None
        if db_key in met.annotations.keys():
            for entry in met.annotations[db_key]:
                new_met_anno_entry, new_names_entry, inchi, new_formula = annotateEntry(entry, db)
                for key, value in new_met_anno_entry.items():
                    if key in new_met_anno.keys():
                        new_met_anno[key].extend(value)
                    else:
                        new_met_anno[key] = value
                # combine the names for each entry
                new_names.extend(new_names_entry)
                # get the optimal inchis
                new_inchi = findOptimalInchi([new_inchi, inchi])

            # add new names to the MeMoMetabolite
            if source is None: raise Exception("source is None")
            if len(new_names) > 0:
                x =met.add_names(new_names, source)
                new_names_added = new_names_added + x
            
            # add the annotations to the slot in the metabolites
            if len(new_met_anno) > 0:
                x = met.add_annotations(new_met_anno, source)
                new_annos_added = new_annos_added + x
            
            # add inchi
            if new_inchi != "":
                x = met.add_inchi_string(new_inchi, source)
                new_inchis_added = new_inchis_added + x

            if new_formula is not None: 
                 x = met.set_formula(new_formula, source)
                 new_formulas_added = new_formulas_added + x

    anno_result = AnnotationResult(new_inchis_added, new_names_added, new_names_added, new_formulas_added)
    return anno_result

