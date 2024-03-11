import urllib.request
import os
from pathlib import Path
from urllib.error import URLError

#from src.config import *
import yaml
import logging

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

def create_folder(path: Path) -> Path:
    """
    param path: Folder where we want to create the new Database folder.
    """
    print("Creating database folder in " + str(path))
    try:
        path.mkdir(parents = True)
    except FileExistsError:
        # if it is an empty directory, that is fine, if not return nothing
        if not len(os.listdir(path)) == 0:
            raise FileExistsError(str(path) +" already exists and is not empty. Please check your database configuration")


def _download(path: Path, URL: str):
    """
    param path: Path to a folder where the db should be stored, in our case this should the /Databases/ folder
    param db: which db to download
    return:
    """
    print("Downloading " + URL + " to " + str(path))
    try:
        urllib.request.urlretrieve(URL, path)
    except URLError:
        print("WARNING: Could not download " + path)


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

    import time
    from concurrent.futures import ThreadPoolExecutor
    start_time =  time.time()

    with ThreadPoolExecutor(max_workers=8) as executor:
      for i in config["databases"]:
          db_path = database_path.joinpath(config["databases"][i]["file"])
          if os.path.exists(db_path):
              os.remove(db_path)
          executor.submit(_download, db_path, config["databases"][i]["URL"])

    # Wait until all downloads are finished
    executor.shutdown(wait=True)
    end_time =  time.time()
    logging.debug(f"Downloading databases took: {end_time - start_time}")

    if "VMH" in config["databases"].keys():
        # handle special case for vmh
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
    return True
