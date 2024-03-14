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
dat = this_directory.joinpath("manually_merged_models")

#class Test_annotateBulkRoutines(unittest.TestCase):
def test_MeMoModelCompare2():
    import time
    
    # Record the start time
    start_time = time.time()
    # test the comparison for metabolite matching
    mod_path = dat.joinpath("MYb11_iCEL1314/M1_MYb11.xml")
    mod_path2 = dat.joinpath("MYb11_iCEL1314/M2_iCEL1314.xml")
    mod = MeMoModel.fromPath(mod_path)
    mod2 = MeMoModel.fromPath(mod_path2)

    print(len(mod.cobra_model.metabolites))
    print(len(mod2.cobra_model.metabolites))

    # self comparison
    res = mod.match(mod2)
    res.to_csv("stuff.csv")
    end_time = time.time()
    
    # Calculate the elapsed time
    elapsed_time = end_time - start_time
    
    print("Elapsed time:", elapsed_time, "seconds")




if __name__ == '__main__':
    test_MeMoModelCompare2()
