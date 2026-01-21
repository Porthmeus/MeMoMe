
#!/usr/bin/env python3

"""
Main entry point of the program
"""
import argparse
import logging
import sys

import cobra

import src.MeMoModel
from src.dbhandling.reformatBiGG import reformatBiGG
from src.dbhandling.reformatModelSeed import reformatModelSeed
from src.dbhandling.reformatHMDB import *
from src.dbhandling.reformatVMH import *
from src.MeMoModel import *
from src.download_db import download, databases_available, update_database
from src.dbhandling.reformatAux import getData, writeData

# Configure the logger
logging.basicConfig(
    level=logging.DEBUG,  # Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# Create a logger with the desired name
logger = logging.getLogger('logger')

# Create a FileHandler to write logs to a file
file_handler = logging.FileHandler('app.log', mode='w')
logger.addHandler(file_handler)


def main(args: argparse.Namespace):
    if args.download:
        logger.debug("Starting to download databases")
        # check if the path database folder exists
        if not databases_available("file"):
          download("URL", "file")
        else:
            update_database("URL", "file")
        logger.debug("Finished downloading databases")
    elif args.reformat:
        dat = getData("BiGG")
        bdb = reformatBiGG(dat)
        writeData(bdb, db = "BiGG")
        
        dat = getData("ModelSeed")
        mdb = reformatModelSeed(dat)
        writeData(mdb, db = "ModelSeed")


        df = getData("HMDB")
        dat_all = do_all(df)
        dat_all.to_csv("hmdb_metabolites_reformatted_all.csv")
        writeData(dat_all, db = "HMDB")


        refDB = reformatVMH()
        writeData(refDB, db = "VMH")


        # Load the model
if __name__ == '__main__':
  # Specifies which arguments are accepted by the program
    parser = argparse.ArgumentParser(description='MeMoMe - Cool stuff.')
    # Specifying this tells the program to download all the databases
    parser.add_argument('--download', action='store_true', help='Download all required databases')
    parser.add_argument('--reformat', action='store_true', help='Reformat all required databases')
    args = parser.parse_args()
    # Log arguments
    logger.debug(args)
    main(args)
