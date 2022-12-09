import sys
import os
import urllib.request
from pathlib import Path
from enum import Enum
from typing import Dict


class Database(Enum):
    ModelSeed = "modelSeed.tsv"
    BiGG = "BiGG.tsv"


database_links: Dict[Database, str] = {
    Database.ModelSeed: "https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/compounds.tsv",
    Database.BiGG: "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt"
}


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


def _download(path: Path, db: Database):
    """
    :param path: Path to a folder where the db should be stored, in our case this should the /Databases/ folder
    :param db: which db to download
    :return:
    """
    print("Downloading + " + str(db.value) + " +  to " + str(path))
    # TODO Find out what happens when retrieving the file fails and handle it
    urllib.request.urlretrieve(database_links[db], path.joinpath(Path(str(db.value))))


def download() -> bool:
    """
    This method creates a folder to store the databases and also downloads them
    :return:
    """
    # Get current path of this script. This should always end in /src and we will create the Database in
    # ../src i.e. in the MeMoMe main directory
    this_path: Path = Path(__file__)
    # path where we want to save the downloaded files
    database_path: Path = create_folder(this_path.parent.parent)
    if database_path == Path(""):
        # Error: We could not create the folder
        return False
    else:
        for k, v in database_links.items():
            _download(database_path, k)
