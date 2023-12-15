#!/usr/bin/env python3

"""
Main entry point of the program
"""
import argparse
import logging
import sys
from os import makedirs os.makedirs
import cobra
from os import makedirs

import src.MeMoModel
from src.MeMoModel import *
from src.download_db import download
# Import the function for duplicate metabolites handling
# from src.removeDuplMetsRxns import update_model_duplicateMetId

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

    # Validate and load both models
    if args.model1 is None or args.model2 is None:
        print("Please supply two models with --model1 and --model2 parameters")
        sys.exit(1)

    model1 = MeMoModel.fromPath(Path(args.model1))
    model2 = MeMoModel.fromPath(Path(args.model2))
    
    v1 = cobra.io.sbml.validate_sbml_model(args.model1)
    v2 = cobra.io.sbml.validate_sbml_model(args.model2)
    print("Validation results for Model 1:", v1)
    print("Validation results for Model 2:", v2)

    # Optional execution of update_model_duplicateMetId based on a feature flag
    if args.update_dup_mets:
        logger.debug("Updating models for duplicate metabolite IDs")
        # update_model_duplicateMetId()  # Uncomment when function is available
        pass


    ann1 = model1.annotate()
    ann2 = model2.annotate()
    
    matches_table = model1.match(model2)
    
    inchi_matches = matches_table[matches_table["inchi_score"]==1]
    makedirs('Output/matches', exist_ok=True)
    inchi_matches.loc[:,["met_id1","met_id2"]].to_csv("Output/matches/inchi_matches.csv")
    unsure_matches = matches_table[matches_table["inchi_score"]!=1]
    unsure_matches.to_csv("Output/matches/to_be_confirmed.csv")

    print(str(len(inchi_matches)) + " metabolite couples have a perfect Inchi string match. The list of matches have been saved under 'Output/matches/inchi_matches.csv\nAdditionally, " +str(len(unsure_matches)) + " potential matches are listed under 'Output/matches/to_be_confirmed.csv'. You can ")
    

if __name__ == '__main__':
    # Specifies which arguments are accepted by the program
    parser = argparse.ArgumentParser(description='MeMoMe - Cool stuff.')
    # Specifying this tells the program to download all the databases
    parser.add_argument('--download', action='store_true', help='Download all required databases')
    parser.add_argument('--model1', action='store', help='Path to the first model that should be merged')
    parser.add_argument('--model2', action='store', help='Path to the second model that should be merged')
    parser.add_argument('--update-dup-mets', action='store_true', help='Update models for duplicate metabolite IDs')
    args = parser.parse_args()
    # Log arguments
    logger.debug(args)
    main(args)
