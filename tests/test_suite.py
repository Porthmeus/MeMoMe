import unittest
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '.')))
print(sys.path)
from matchInchi import *
from annotate_dbs import *
from annotateInchiRoutines import *
from matchInchi import *
from MeMoMetabolite import *
from MeMoModel import *
from parseMetaboliteInfoFromSBML import *
from PubchemInfoNoInet import *
from PubchemInfo import *
from removeDuplicateMetabolites import *
from bulkPerformance import *



fast_tests = [Test_annotateMissingDbs,Test_annotateEntryFunctions,Test_annotateID,Test_annotateFull,TestInchiRoutines,TestMatchInchit,Test_MeMoMetabolite,Test_annotateBulkRoutines,Test_MiscStuff,Test_removeDuplicates,TestParseSBML,TestPubchemInfoNoInet,TestPubchemInfo,Test_removeDuplicateMetabolites]


slow_tests = [Test_annotationPerformance]

def fast():
  suite = unittest.TestSuite()
  for x in fast_tests:
    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(x))
  return suite

def slow():
  suite = unittest.TestSuite()
  for x in slow_tests:
    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(x))
  return suite

if __name__ == '__main__':
  runner = unittest.TextTestRunner()
  if len(sys.argv) < 2:
    print("Please specify fast or slow")

  fasts = fast()
  slows = slow()

  if sys.argv[1] == "fast":
    print("\nRunning Unit Tests...")
    runner.run(fasts)

  if sys.argv[1] == "slow":
    print("\nRunning System Tests...")
    runner.run(slows)
