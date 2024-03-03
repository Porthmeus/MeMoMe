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
from src.download_db import download

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
        download()
        logger.debug("Finished downloading databases")

    # Check if exactly two models were supplied
    if args.model1 is None:
        print("Please supply a first model with the --model1 parameter", file=sys.stderr)
        logger.error("User did not provide the path to the first model")
        sys.exit(1)
    if args.model2 is None:
        print("Please supply a second model with the --model2 parameter", file=sys.stderr)
        logger.error("User did not provide the path to the second model")
        sys.exit(1)
    if args.output is None:
        print("Please provide at path and output file name <path>/<outname>.csv", file=sys.stderr)
        logger.error("User did not provide an output path")
        sys.exit(1)
    
    # Load the model
    model1 = MeMoModel.fromPath(Path(args.model1))
    model2 = MeMoModel.fromPath(Path(args.model2))
    # bulk annotate the model
    model1.annotate()
    model2.annotate()

    matched_model = model1.match(model2)

    matched_model.to_csv(args.output, index = False)

if __name__ == '__main__':
    # Specifies which arguments are accepted by the program
    parser = argparse.ArgumentParser(description='MeMoMe - Cool stuff.')
    # Specifying this tells the program to download all the databases
    parser.add_argument('--download', action='store_true', help='Download all required databases')
    parser.add_argument('--model1', action='store', help='Path to the first model that should be merged')
    parser.add_argument('--model2', action='store', help='Path to the second model that should be merged')
    parser.add_argument('--output', action='store', help='Path where the output should be stored (as a csv)')
    args = parser.parse_args()
    # Log arguments
    logger.debug(args)
    main(args)
