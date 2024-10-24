# Porthmeus
# 08.03.24

from __future__ import annotations
#### original code by @unaimed ####
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
