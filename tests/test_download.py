import unittest
from pathlib import Path

from src import download_db


class TestDownload(unittest.TestCase):
    # The directory of this file
    this_directory = Path(__file__).parent

    def test_non_ex_parent(self):
        # create a nonexistent target to download the stuff to
        non_exist_path = self.this_directory.parent.joinpath(Path("nonexist"))
        # The call to create_folder should return Path("")
        self.assertEqual(download_db.create_folder(non_exist_path), Path(""))

    def test_dir_alr_exists(self):
        self.assertEqual(download_db.create_folder(self.this_directory), self.this_directory.joinpath("Databases"))
        self.assertEqual(download_db.create_folder(self.this_directory), Path(""))

        # Because other methods might try to create this folder we will delete it again
        self.this_directory.joinpath("Databases").rmdir()

    def test_exist_parent(self):
        # Check that create folder creates a /Database/ folder
        self.assertEqual(download_db.create_folder(self.this_directory), self.this_directory.joinpath("Databases"))

    @classmethod
    def tearDownClass(cls):
        # Delete the created folder for the test.
        # This is done after all tests have been executed.
        cls.this_directory.joinpath("Databases").rmdir()


if __name__ == '__main__':
    unittest.main()
