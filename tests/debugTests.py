import unittest

class Test_debug(unittest.TestCase):
    def test_debug_fail():
        raise Exception("Debug fail")
