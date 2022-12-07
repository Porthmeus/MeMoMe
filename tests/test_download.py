import unittest
from pathlib import Path

from src import download_db


class TestDownload(unittest.TestCase):
    def test_non_ex_parent(self):
        # get current path
        this_path = Path(__file__)
        # create a nonexistent target to download the stuff to
        non_exist_path = this_path.parent.joinpath(Path("nonexist"))
        # The call to creat_folder should return Path("")
        self.assertEqual(download_db.create_folder(non_exist_path), Path(""))


if __name__ == '__main__':
    unittest.main()
