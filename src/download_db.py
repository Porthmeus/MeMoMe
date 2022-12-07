import sys
import os
import urllib.request
from pathlib import Path

def create_folder(path: Path) -> Path:
    """
    :param path: Folder where we want to create the new Database folder
    :return: "" if no folder was created, otherwise the path to the folder
    """
    print("Creating database folder in " + str(path))
    db_path = path.joinpath("Databases")
    # TODO: Handle file exists
    try:
        db_path.mkdir()
    except FileExistsError:
        return Path("")
    except FileNotFoundError:
        return Path("")

    return db_path


def download_modelseed(path: Path):
    # TODO Findout what happens when retrieving the file fails and handle it
    urllib.request.urlretrieve(
        "https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/compounds.tsv",
        path.joinpath(Path("modelSeed.tsv")))


def download() -> bool:
    """
    This method creates a folder to store the databases and also downloads them
    :return:
    """
    # Get current path of this script. This should always end in /src and we will create the Database in
    # ../src i.e. in the MeMoMe main directory
    this_path: Path = Path(__file__)
    # path where we want to save the downloaded files
    database_path = create_folder(this_path.parent.parent)
    if database_path == "":
        # Error: We could not create the folder
        return False
    else:
        download_modelseed(database_path)
