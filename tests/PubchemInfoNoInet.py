# Porthmeus
# 14.05.24

# block internet connection for python
import socket


class block_network(socket.socket):
    def __init__(self, *args, **kwargs):
        raise socket.timeout("Network call blocked")


socket.socket = block_network


import unittest
from pathlib import Path
import pandas as pd
from src.downloadPubchemInfo import looksLikeSumFormula, getPubchemInfo


class TestPubchemInfoNoInet(unittest.TestCase):
    # The directory of this file
    this_directory = Path(
        __file__
    ).parent  # parent because the actual file is part of the path
    dat_dir = this_directory.joinpath("dat")

    def test_getPubchemInfo_noInternet(self):
        self.assertIsInstance(getPubchemInfo(["alkdjfa", "Carnitine"]), pd.DataFrame)
        self.assertTrue(getPubchemInfo(["alkdjfa", "Carnitine"]).empty)


if __name__ == "__main__":
    unittest.main()
