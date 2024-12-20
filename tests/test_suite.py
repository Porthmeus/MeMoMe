import unittest
from test_matchInchi import *
from test_annotate_dbs import *
from test_annotateInchiRoutines import *
from test_matchInchi import *
from test_MeMoMetabolite import *
from test_MeMoModel import *
from test_parseMetaboliteInfoFromSBML import *
from test_PubchemInfoNoInet import *
from test_PubchemInfo import *
from test_removeDuplicateMetabolites import *



fast_tests = [Test_annotateMissingDbs,Test_annotateEntryFunctions,Test_annotateID,Test_annotateFull,TestInchiRoutines,TestMatchInchit,Test_MeMoMetabolite,Test_annotateBulkRoutines,Test_MiscStuff,Test_removeDuplicates,TestParseSBML,TestPubchemInfoNoInet,TestPubchemInfo,Test_removeDuplicateMetabolites]


def fast():
  suite = unittest.TestSuite()
  for x in fast_tests:
    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(x))
  return suite


#def slow():
#  suite = unittest.TestSuite()
#  suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestSystemLevel))
#  return suite

if __name__ == '__main__':
  runner = unittest.TextTestRunner()

  print("\nRunning Unit Tests...")
  runner.run(fast())

  #print("\nRunning System Tests...")
  #runner.run(system_suite())
