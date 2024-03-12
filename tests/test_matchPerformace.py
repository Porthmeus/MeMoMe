# Porthmeus
# 02.08.23

import unittest
from pathlib import Path
import cobra as cb
import pandas as pd
import warnings

# from src.MeMoModel import MeMoModel
from src.MeMoMetabolite import MeMoMetabolite
from src.MeMoModel import MeMoModel

this_directory = Path(__file__).parent
dat = this_directory.joinpath("dat")

#class Test_annotateBulkRoutines(unittest.TestCase):
def test_MeMoModelCompare2():
    import time
    
    # Record the start time
    start_time = time.time()
    # test the comparison for metabolite matching
    mod_path = dat.joinpath("Recon3DModel_301.xml")
    mod = MeMoModel.fromPath(mod_path)
    # self comparison
    res = mod.match(mod)

    end_time = time.time()
    
    # Calculate the elapsed time
    elapsed_time = end_time - start_time
    
    print("Elapsed time:", elapsed_time, "seconds")




if __name__ == '__main__':
    test_MeMoModelCompare2()
