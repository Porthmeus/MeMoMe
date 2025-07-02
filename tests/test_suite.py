import unittest
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

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
from DBreformatting import *


fast_tests = [Test_annotateMissingDbs,Test_annotateEntryFunctions,Test_annotateID,Test_annotateFull,TestInchiRoutines,TestMatchInchit,Test_MeMoMetabolite,Test_annotateBulkRoutines,Test_MiscStuff,Test_removeDuplicates,TestParseSBML,TestPubchemInfoNoInet,TestPubchemInfo,Test_removeDuplicateMetabolites]


slow_tests = [Test_annotationPerformance]

formatting_tests = [TestConcatCols, TestRenameColumnsSafe, TestPrepare, TestGetAnnos, TestChebi]

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

def formatting():
  suite = unittest.TestSuite()
  for x in formatting_tests:
    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(x))
  return suite

if __name__ == '__main__':
  runner = unittest.TextTestRunner(stream=sys.stdout, verbosity=2, buffer=False)
  if len(sys.argv) < 2:
    print("Please specify fast or slow")

  fasts = fast()
  slows = slow()
  formattings = formatting()
  if sys.argv[1] == "fast":
    print("\nRunning Unit Tests...")
    runner.run(fasts)

  if sys.argv[1] == "slow":
    print("\nRunning System Tests...")
    runner.run(slows)

  if sys.argv[1] == "formatting":
    print("\nRunning Formmating Tests...")
    runner.run(
        formattings
        )
