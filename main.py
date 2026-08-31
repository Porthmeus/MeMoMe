#!/usr/bin/env python3

"""
Main entry point of the programt
"""
import argparse
import logging
import sys
import cobra
import os

from pathlib import Path
from src.MeMoModel import *
from src.ModelMerger import ModelMerger
from src.download_db import download, databases_available, update_database
from src.matchMets import filter_matching_table


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
        if not databases_available("reformat"):
            download("REFORMAT_URL", "reformat")
        else:
            update_database("REFORMAT_URL", "reformat")
        logger.debug("Finished downloading databases")
    else:
        # Check if exactly two models were supplied
        if args.model1 is None:
            print("Please supply a second model with the --model1 parameter")
            sys.exit(1)
        if args.model2 is None:
            print("Please supply a second model with the --model2 parameter")
            sys.exit(1)


        # Load the model
        model1 = MeMoModel.fromPath(Path(args.model1))
        model2 = MeMoModel.fromPath(Path(args.model2))
        # bulk annotate the model
        model1.annotate(args.allow_missing_dbs)
        model2.annotate(args.allow_missing_dbs)
        # create output directory
        output_path = args.output
        output_path.mkdir(parents=True, exist_ok=True)
        # do the matching and safe the matching table
        matching_table = model1.match(model2, output_names = args.output_names, output_dbs = args.output_dbs, keepUnmatched = args.keep_unmatched)
        matching_table.to_csv(args.output/Path("matching_table.csv"), index = False)
        # filter matching table
        matching_table = filter_matching_table(matching_table,
                                               Inchi_threshold = args.InChI_threshold,
                                               DB_threshold = args.DB_threshold,
                                               Name_threshold = args.Name_threshold)

        # merge the models and save to sbml
        merger = ModelMerger(model1, model2, matching_table)
        merger.preprocess_models()
        merger.translate()
        merger.merge_models()
        print(type(args.merged_output))
        cobra.io.write_sbml_model(merger.merged_model,
                                  output_path / args.merged_output)
        if args.save_translated_model:
            split_models = merger.split_merged_model()
            for i, nm  in enumerate(merger.merged_model.notes["MeMoMe_prefixes"].keys()):
                cobra.io.write_sbml_model(split_models[i],
                                          output_path / Path(nm+".sbml"))



if __name__ == '__main__':
    # Specifies which arguments are accepted by the program
    parser = argparse.ArgumentParser(description='MeMoMe - Cool stuff.')
    # Specifying this tells the program to download all the databases
    parser.add_argument('--output-names', action='store_true', default=False, help='If two metabolites got matched on a name basis, output the names that lead to this match')
    parser.add_argument('--output-dbs', action='store_true', default=False, help='If two metabolites got matched on a database basis, output the databases that lead to this match')
    parser.add_argument('--keep-unmatched', action='store_true', default=False, help='Stored unmatched metabolties in the output')
    parser.add_argument('--download', action='store_true', help='Download all required databases')
    parser.add_argument('--reformat', action='store_true', help='Reformat all required databases')
    parser.add_argument('--model1', action='store', help='Path to the first model that should be merged')
    parser.add_argument('--model2', action='store', help='Path to the second model that should be merged')
    parser.add_argument('--output', action='store', default = "MeMoMe_output", help='Path where of the output directory', type = Path)
    parser.add_argument('--allow_missing_dbs', action='store_true', help='If set to true program does not abort if a databse is missing')
    parser.add_argument('--merged-output', action='store', default = "merged_model.xml", help='Path where the merged model should be stored (as an SBML file)', type = Path)
    parser.add_argument('--save-translated-model', action='store_true', default=False, help='Should the merged models be saved as individual sbml files?')
    parser.add_argument("--InChI-threshold", action = "store", default = 0, type = float, help = "Minimum InChI threshold which needs to be achieved to retain a match in the matching table (note InChI score is either 0 or 1)")
    parser.add_argument("--DB-threshold", action = "store", default = 0, type = float, help = "Minimum DB threshold which needs to be achieved to retain a match in the matching table")
    parser.add_argument("--Name-threshold", action = "store", default = 0, type = float, help = "Minimum name threshold which needs to be achieved to retain a match in the matching table")

    args = parser.parse_args()
    # Log arguments
    logger.debug(args)
    main(args)
