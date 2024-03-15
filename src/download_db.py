import urllib.request
import os
from pathlib import Path
from urllib.error import URLError

#from src.config import *
import yaml
import requests
import sys
def is_inet_available() -> bool:
    try:
        response = requests.get("http://www.google.com", timeout=5)
        response.raise_for_status()  # Raises an exception for 4xx and 5xx status codes
        print("Internet is available.")
        return True
    except requests.RequestException:
        print("Internet is not available.")
        return False

def get_config() -> dict:
    """ load the config file and return the dictionary"""
    this_path: Path = Path(__file__)
    config_path: Path = this_path.parent.parent.joinpath("config.yml")
    with open(config_path, "r") as ff:
        config = yaml.safe_load(ff)
    return config


def get_database_path() -> Path:
    """
    Simple wrapper to extract the database path from the config file. Returns an absolute Path to the database location
    """
    this_path: Path = Path(__file__)
    config_path: Path = this_path.parent.parent.joinpath("config.yml")
    with open(config_path, "r") as ff:
        config = yaml.safe_load(ff)
    database_path:Path = this_path.parent.parent.joinpath(config["DB_location"])
    return(database_path.absolute())

def create_folder(path: Path) -> Path | None:
    """
    param path: Folder where we want to create the new Database folder.
    """
    try:
        dbfolder_there  = path.exists()
        inet_there = is_inet_available()
        if (not dbfolder_there ) and (not inet_there):
          print("No databases were found but no internet connection is available, aborting program")
          sys.exit(1)
        # TODO This code is not testable
        print("Creating database folder in " + str(path))
        path.mkdir(parents = True)
    except FileExistsError:
        # if it is an empty directory, that is fine, if not return nothing
        if not len(os.listdir(path)) == 0:
            raise FileExistsError(str(path) + f" already exists and is not empty. Please check your database configuration or delete content in {str(path)})")


def _download(path: Path, URL: str):
    """
    param path: Path to a folder where the db should be stored, in our case this should the /Databases/ folder
    param db: which db to download
    return:
    """
    try:
        urllib.request.urlretrieve(URL, path)
        print("Downloading " + URL + " to " + str(path))
    except URLError:
        print("WARNING: Could not download " + str(path))


def download() -> bool:
    """
    This method creates a folder to store the databases and also downloads them
    :return:
    """
    # Get current path of this script. This should always end in /src and we will create the Database as defined in the config file 
    # load the config yaml
    config = get_config()
    database_path = get_database_path()
    # path where we want to save the downloaded files
    create_folder(database_path)

    status = update_database()

    return status

def databases_available() -> bool:
    """
    Method just checks if database folder is present and if all databases have been downloaded - return True if everything is in place, False, if something is missing
    """
    is_there = True

    config = get_config()
    database_path = get_database_path()
    # check if all database files exists if not try to download
    if not os.path.exists(database_path):
        is_there = download()
    else:
        for db in config["databases"]:
            if not os.path.exists(Path(database_path, config["databases"][db]["file"])):
                is_there = False
                break
    return is_there

def update_database() -> bool:
    '''Remove the databases and download them again with possible new updates'''
    # load the config yaml
    config = get_config()
    database_path = get_database_path()

    for i in config["databases"]:
        db_path = database_path.joinpath(config["databases"][i]["file"])
        if os.path.exists(db_path):
            os.remove(db_path)
        _download(db_path, config["databases"][i]["URL"])
    if "VMH" in config["databases"].keys():
        # handle special case for vmh
        try:
          with open(database_path.joinpath(config["databases"]["VMH"]["file"]), mode='r+') as f:
            content: str = f.read()
            # This will break if the link changes
            content = content.removeprefix("Ext.data.JsonP.callback19(")
            content = content.removesuffix(");")
            # Go to the beginning of the file
            f.seek(0)
            # Write the new content
            f.write(content)
            # Truncate to the new contents length(because the old content of the file was longer)
            f.truncate()
        except FileNotFoundError: 
          print("VMH could not be found, not tryting to reformat it")
    return True
