#!/usr/bin/env python3

"""
Main entry point of the program
"""
import argparse
import logging
import sys

import cobra

import src.MeMoModel
from src.MeMoModel import *
from src.download_db import download, databases_available, update_database

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
        if not databases_available():
            download()
        else:
            update_database()
        logger.debug("Finished downloading databases")
    else:
        # Check if exactly two models were supplied
        if args.model1 is None:
            print("Please supply a second model with the --model1 parameter")
            sys.exit(1)
        if args.model2 is None:
            print("Please supply a second model with the --model2 parameter")
            sys.exit(1)

        v = cobra.io.sbml.validate_sbml_model(args.model2)
        print(v)
        # model1 = MeMoModel.fromPath(Path(args.model1))
        model2 = MeMoModel.fromPath(Path(args.model2))

        #t = model1.annotate()
        t = model2.annotate()
        print("T")




if __name__ == '__main__':
    # Specifies which arguments are accepted by the program
    parser = argparse.ArgumentParser(description='MeMoMe - Cool stuff.')
    # Specifying this tells the program to download all the databases
    parser.add_argument('--download', action='store_true', help='Download all required databases')
    parser.add_argument('--model1', action='store', help='Path to the first model that should be merged')
    parser.add_argument('--model2', action='store', help='Path to the second model that should be merged')
    args = parser.parse_args()
    # Log arguments
    logger.debug(args)
    main(args)
