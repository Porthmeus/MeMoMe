# Porthmeus
# 24.02.23

import unittest
from pathlib import Path
import pandas as pd

from src.downloadPubchemInfo import looksLikeSumFormula, getPubchemInfo


class TestPubchemInfo(unittest.TestCase):
    # The directory of this file
    this_directory = Path(__file__).parent # parent because the actual file is part of the path
    dat_dir = this_directory.joinpath("dat")

    def test_looksLikeSumFormula(self):
        self.assertTrue(looksLikeSumFormula("C7H15NO3"))
        self.assertFalse(looksLikeSumFormula("Carnitine"))

    def test_getPubchemInfo(self):
        self.assertIsInstance(getPubchemInfo(["alkdjfa", "Carnitine"]), pd.DataFrame)

    def test_getPubchemInfo_noInternet(self):
        self.assertIsInstance(getPubchemInfo(["alkdjfa", "Carnitine"]), pd.DataFrame)
        self.assertTrue(getPubchemInfo(["alkdjfa", "Carnitine"]).empty)

if __name__ == "__main__":
    unittest.main()
