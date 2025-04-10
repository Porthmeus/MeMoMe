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
from src.TranslationAdder import ModelMerger
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
        # Check if exactly two models were supplied
        if args.output is None:
            print("Please provide at path and output file name <path>/<outname>.csv", file=sys.stderr)
            logger.error("User did not provide an output path")
            sys.exit(1)
    
        # Load the model
        model1 = MeMoModel.fromPath(Path(args.model1))
        model2 = MeMoModel.fromPath(Path(args.model2))
        # bulk annotate the model
        model1.annotate(args.allow_missing_dbs)
        model2.annotate(args.allow_missing_dbs)

        matched_model = model1.match(model2, keep1ToMany = args.keep_one_to_many, output_names = args.output_names, output_dbs = args.output_dbs, keepUnmatched = args.keep_unmatched)
        matched_model.to_csv(args.output, index = False)
        merger = ModelMerger(model2, matched_model)
        # add translation compartment to model2, so that its namespace fits the namespace of model1
        merger.translate_namespace()
        cobra.io.write_sbml_model(model2.cobra_model, args.output_translated_model)

if __name__ == '__main__':
    # Specifies which arguments are accepted by the program
    parser = argparse.ArgumentParser(description='MeMoMe - Cool stuff.')
    # Specifying this tells the program to download all the databases
    parser.add_argument('--output-names', action='store_true', default=False, help='If two metabolites got matched on a name basis, output the names that lead to this match')
    parser.add_argument('--output-dbs', action='store_true', default=False, help='If two metabolites got matched on a database basis, output the databases that lead to this match')
    parser.add_argument('--keep-one-to-many', action='store_true', default=False, help='Keep one-to-many merges')
    parser.add_argument('--keep-unmatched', action='store_true', default=False, help='Stored unmatched metabolties in the output')
    parser.add_argument('--download', action='store_true', help='Download all required databases')
    parser.add_argument('--model1', action='store', help='Path to the first model that should be merged')
    parser.add_argument('--model2', action='store', help='Path to the second model that should be merged')
    parser.add_argument('--output', action='store', help='Path where the output should be stored (as a csv)')
    parser.add_argument('--output_translated_model', action='store', help='Path where model2 with the newly-added translation compartment should be stored (as an SBML file)')
    parser.add_argument('--allow_missing_dbs', action='store_true', help='If set to true program does not abort if a databse is missing')
    args = parser.parse_args()
    # Log arguments
    logger.debug(args)
    main(args)
