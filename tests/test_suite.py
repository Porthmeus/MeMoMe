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
from ModelMerger import *
from bulkPerformance import *
from DBreformatting import *
from debugTests import *


fast_tests = [Test_annotateMissingDbs,
              Test_annotateEntryFunctions,
              Test_annotateID,
              Test_annotateFull,
              TestInchiRoutines,
              TestMatchInchit,
              Test_MeMoMetabolite,
              Test_annotateBulkRoutines,
              Test_MiscStuff,
              Test_ModelMerger,
              Test_removeDuplicates,
              TestParseSBML,
              TestPubchemInfoNoInet,
              TestPubchemInfo,
              Test_removeDuplicateMetabolites
              ]


slow_tests = [
  Test_annotationPerformance, 
  Test_ModelMerger_Slow
  ]

formatting_tests = [TestConcatCols, TestRenameColumnsSafe, TestPrepare, TestGetAnnos, TestChebi]

debug_tests = [Test_debug]

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

def debug():
  suite = unittest.TestSuite()
  for x in debug_tests:
    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(x))
  return suite

if __name__ == '__main__':
  runner = unittest.TextTestRunner(stream=sys.stdout, verbosity=2, buffer=False)
  if len(sys.argv) < 2:
    print("Please specify fast, slow, or formatting")

  fasts = fast()
  slows = slow()
  formattings = formatting()
  debug = debug()
  if sys.argv[1] == "fast":
    print("\nRunning Unit Tests...")
    tests = runner.run(fasts)
    if len(tests.errors) > 0 or len(tests.failures) > 0:
        raise Exception("tests:{run} errors:{errors} failures:{failures}".format(run = tests.testsRun, errors = len(tests.errors), failures = len(tests.failures)))

  if sys.argv[1] == "slow":
    print("\nRunning System Tests...")
    tests = runner.run(slows)
    if len(tests.errors) > 0 or len(tests.failures) > 0:
        raise Exception("tests:{run} errors:{errors} failures:{failures}".format(run = tests.testsRun, errors = len(tests.errors), failures = len(tests.failures)))

  if sys.argv[1] == "formatting":
    print("\nRunning Formmating Tests...")
    tests = runner.run(formattings)
    if len(tests.errors) > 0 or len(tests.failures) > 0:
        raise Exception("tests:{run} errors:{errors} failures:{failures}".format(run = tests.testsRun, errors = len(tests.errors), failures = len(tests.failures)))

  if sys.argv[1] == "debug":
    print("\nRunning debugging Tests...")
    tests = runner.run(debug)
    if len(tests.errors) > 0 or len(tests.failures) > 0:
        raise Exception("tests:{run} errors:{errors} failures:{failures}".format(run = tests.testsRun, errors = len(tests.errors), failures = len(tests.failures)))
