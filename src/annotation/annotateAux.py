# Porthmeus
# 08.03.24
from __future__ import annotations
import pandas as pd
from src.download_db import get_config, get_database_path
import warnings
import os
from typing import Callable, List
from src.MeMoMetabolite import MeMoMetabolite

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
                  conversion_method: Callable[[str], pd.DataFrame] = lambda x: x) -> pd.DataFrame:
  """
  Load the given database. The file should in the projects root /Database folder.
  """
  # lad the database
  try:
    db_path =  os.path.join(get_database_path(), database)
    db = conversion_method(db_path)
  except FileNotFoundError as e:
    warnings.warn(str(e))
    # Rethrow exception because we want don't allow missing dbs
    if allow_missing_dbs == False:
      raise e
    db = pd.DataFrame()
  return(db)

def handleIDs(db: pd.DataFrame, metabolites: List[MeMoMetabolite], db_key: str, annotation_function: Callable[[str, pd.DataFrame], tuple[dict, list]]) -> AnnotationResult:
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
      new_met_anno_entry,new_names_entry = annotation_function(met._id, db)
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


  
