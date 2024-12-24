# Porthmeus
# 08.03.24
from __future__ import annotations
import pandas as pd
from src.download_db import get_config, get_database_path
import warnings
import os
from typing import Callable, List, NewType
from src.MeMoMetabolite import MeMoMetabolite
import sys

# HMDB0000972
AnnotationKey = NewType('AnnotationKey', str)
# Database name, corresponds to one of the keys in config.yaml
DBName = NewType('DBName', str)
# The dbs have have certain columns that correspond to their main identifier, e.g. vmh has vmhmetabolite, HMDB has accession and so one. This is what we call a database key. Aka the MAIN column which we use for entry lookup
DBKey = NewType('DBKey', str)

EntryAnnotationFunction = Callable[[AnnotationKey, pd.DataFrame, bool], tuple[dict, list]]



class AnnotationResult():
  def __init__(self, annotated_inchis: int, annotated_dbs: int, annotated_names: int):
    self.annotated_inchis: int =  annotated_inchis
    self.annotated_dbs: int = annotated_dbs
    self.annotated_names: int = annotated_names
    self.annotated_total: int = annotated_inchis + annotated_dbs + annotated_names
  
  @classmethod
  def fromTuple(cls, annotationResult: tuple[int, int, int]) -> AnnotationResult:
    return cls(annotationResult[0], annotationResult[1], annotationResult[2])
  
  @classmethod
  def fromAnnotation(cls, annotationResult: AnnotationResult) -> AnnotationResult:
    return cls(annotationResult.annotated_inchis, annotationResult.annotated_dbs, annotationResult.annotated_names)

  def __str__(self) -> str:
    return f"Annotated inchis {self.annotated_inchis}, annotated dbs {self.annotated_dbs}, annotated names {self.annotated_names}"

  def __repr__(self) -> str:
    return f"Annotated inchis {self.annotated_inchis}, annotated dbs {self.annotated_dbs}, annotated names {self.annotated_names}"

  def __le__(self, other) -> bool:
    return self.annotated_inchis <= other.annotated_inchis and self.annotated_dbs <= other.annotated_dbs and self.annotated_names <= other.annotated_names

  def __gt__(self, other):
    return not (self.__le__(other))

  def __add__(self, other):
    return AnnotationResult(self.annotated_inchis + other.annotated_inchis, self.annotated_dbs + other.annotated_dbs, self.annotated_names + other.annotated_names)

  def __iter__(self):
      return iter([self.annotated_inchis, self.annotated_dbs, self.annotated_names])

  def __eq__(self, other): 
      if not isinstance(other, AnnotationResult):
          # don't attempt to compare against unrelated types
          return NotImplemented

      return self.annotated_inchis == other.annotated_inchis and self.annotated_dbs == other.annotated_dbs and self.annotated_names == other.annotated_names
####################################


def load_database(database: str = "", allow_missing_dbs: bool = False, 
                  conversion_method: Callable[[str], pd.DataFrame] = lambda x: pd.DataFrame()) -> pd.DataFrame:
  """
  Load the given database. The file should in the projects root /Database folder.
  """
  # load the database
  try:
    db_path =  os.path.join(get_database_path(), database)
    db = conversion_method(db_path)
  except FileNotFoundError as e:
    # Rethrow exception if we don't allow missing dbs
    if allow_missing_dbs == False:
      raise e
    warnings.warn(str(e))
    db = pd.DataFrame()
  return(db)

def handleIDs(db: pd.DataFrame, metabolites: List[MeMoMetabolite], db_key: str, entry_annotation_function: EntryAnnotationFunction, allow_missing_dbs: bool = False) -> AnnotationResult:
  """
  Checks for each metabolite if the metabolite id can be found in the column `db_key` can be found in `db`. 
  db: a dataframe with columns `db_key`. This column will be comparted to the metabolite id (met._id)
  metabolites: A list of metabolites that will be checked
  db_key: The column in the db dataframe
  annotation_function: Defines how to get an entry from the db which the given met._id (Check annotateVMH/BiGG for example usages.
  """
  new_annos = 0
  new_names = 0
  for met in metabolites:
    if any(db[db_key]==met._id):
      #TODO FIX LAST PARaM
      new_met_anno_entry,new_names_entry = entry_annotation_function(met._id, db, allow_missing_dbs)
      x = met.add_names(new_names_entry)
      new_names = new_names + x

      # add the annotations to the slot in the metabolites
      if len(new_met_anno_entry) > 0:
          x = met.add_annotations(new_met_anno_entry)
          new_annos = new_annos + x 

  # return the number of metabolites which got newly annotated with inchis,
  # annotations and names
  anno_result = AnnotationResult(0, new_annos, new_names)
  return anno_result


def handleMetabolites(db: pd.DataFrame, metabolites: List[MeMoMetabolite], db_key: str, annotation_function: EntryAnnotationFunction, allow_missing_dbs: bool = False) -> AnnotationResult:
    new_annos_added = 0
    new_names_added = 0
    # go through the metabolites and check if there is data whigh can be added
    for met in metabolites:
        new_met_anno = dict()
        new_names = list() 
        if db_key in met.annotations.keys():
            for entry in met.annotations[db_key]:
                new_met_anno_entry, new_names_entry = annotation_function(entry, db, allow_missing_dbs)
                for key, value in new_met_anno_entry.items():
                    if key in new_met_anno.keys():
                        new_met_anno[key].extend(value)
                    else:
                        new_met_anno[key] = value
                # combine the names for each entry
                new_names.extend(new_names_entry)

            # add new names to the MeMoMetabolite
            x =met.add_names(new_names)
            new_names_added = new_names_added + x
            
            # add the annotations to the slot in the metabolites
            if len(new_met_anno) > 0:
                x = met.add_annotations(new_met_anno)
                new_annos_added = new_annos_added + x
    anno_result = AnnotationResult(0, new_names_added, new_names_added)
    return anno_result
