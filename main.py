#!/usr/bin/env python3

"""
Main entry point of the programt
"""
import argparse
import logging
import sys
import cobra
import os
import pickle
import pandas as pd

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


def annotate_models(args: argparse.Namespace):
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
    for model, filename in [(model1, "model1.pkl"), (model2, "model2.pkl")]:
        with (args.output / filename).open("wb") as handle:
            pickle.dump(model, handle)


def match_models(args: argparse.Namespace):
    with (args.output / "model1.pkl").open("rb") as handle:
        model1 = pickle.load(handle)
    with (args.output / "model2.pkl").open("rb") as handle:
        model2 = pickle.load(handle)
    # do the matching and safe the matching table
    matching_table = model1.match(model2, output_names = args.output_names, output_dbs = args.output_dbs, keepUnmatched = args.keep_unmatched)
    matching_table.to_csv(args.output / "matching_table.csv", index = False)


def merge_models(args: argparse.Namespace):
    with (args.output / "model1.pkl").open("rb") as handle:
        model1 = pickle.load(handle)
    with (args.output / "model2.pkl").open("rb") as handle:
        model2 = pickle.load(handle)
    matching_table_path = args.matching_table if args.matching_table is not None else args.output / "matching_table.csv"
    matching_table = pd.read_csv(matching_table_path, dtype={"met_id1": str, "met_id2": str})
    # filter matching table
    matching_table = filter_matching_table(matching_table,
                                           Inchi_threshold = args.InChI_threshold,
                                           DB_threshold = args.DB_threshold,
                                           Name_threshold = args.Name_threshold,
                                           score_type = args.score_type)

    # merge the models and save to sbml
    merger = ModelMerger(model1, model2, matching_table)
    merger.preprocess_models()
    merger.translate()
    merger.merge_models()
    cobra.io.write_sbml_model(merger.merged_model,
                              args.output / args.merged_output)
    if args.save_translated_model:
        split_models = merger.split_merged_model()
        for i, nm in enumerate(merger.merged_model.notes["MeMoMe_prefixes"].keys()):
            cobra.io.write_sbml_model(split_models[i],
                                      args.output / Path(nm+".sbml"))


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
        # create output directory
        args.output.mkdir(parents=True, exist_ok=True)
        if args.step in ("all", "annotate"):
            annotate_models(args)
        if args.step in ("all", "match") and not (args.step == "all" and args.matching_table is not None):
            match_models(args)
        if args.step in ("all", "merge"):
            merge_models(args)



if __name__ == '__main__':
    # Specifies which arguments are accepted by the program
    parser = argparse.ArgumentParser(description='MeMoMe - Cool stuff.')
    parser.add_argument('--step', choices=['all', 'annotate', 'match', 'merge'], default='all', help='''Workflow step to run (default: all);
        annotate - only annotates the two models (requires <model1> and <model2>);
        match - matches the metabolites of the two annotated models in <output>;
        merge - merges the two annotated models in <output> with the table provided by <matching_table>, if <matching_table> is not given, default matching table in <output> is used.''')
    parser.add_argument('--matching-table', type=Path, help='CSV matching table to use for merging; skips matching when running all steps (default: creates and uses matching_table.csv in the output directory)')
    # Specifying this tells the program to download all the databases
    parser.add_argument('--output-names', action='store_true', default=False, help='If two metabolites got matched on a name basis, output the names that lead to this match')
    parser.add_argument('--output-dbs', action='store_true', default=False, help='If two metabolites got matched on a database basis, output the databases that lead to this match')
    parser.add_argument('--keep-unmatched', action='store_true', default=False, help='Stored unmatched metabolties in the output')
    parser.add_argument('--download', action='store_true', help='Download all required databases')
    parser.add_argument('--reformat', action='store_true', help='Reformat all required databases')
    parser.add_argument('--model1', action='store', help='Path to the first model that should be merged (SBML format)')
    parser.add_argument('--model2', action='store', help='Path to the second model that should be merged (SBML format)')
    parser.add_argument('--output', action='store', default = "MeMoMe_output", help='Path where of the output directory', type = Path)
    parser.add_argument('--allow_missing_dbs', action='store_true', help='If set to true program does not abort if a databse is missing')
    parser.add_argument('--merged-output', action='store', default = "merged_model.xml", help='Path where the merged model should be stored (as an SBML file)', type = Path)
    parser.add_argument('--save-translated-model', action='store_true', default=False, help='Should the merged models be saved as individual sbml files?')
    parser.add_argument("--InChI-threshold", action = "store", default = 1.0, type = float, help = "Minimum InChI score sufficient to retain a match on its own (default: 1.0; score is either 0 or 1)")
    parser.add_argument("--DB-threshold", action = "store", default = 0.5, type = float, help = "Minimum DB score required together with the name threshold when the InChI threshold is not met (default: 0.5)")
    parser.add_argument("--Name-threshold", action = "store", default = 0.9, type = float, help = "Minimum name score required together with the DB threshold when the InChI threshold is not met (default: 0.9)")

    parser.add_argument("--score-type", default="total_score", help="Score column used to select the best match per metabolite (default: total_score)")

    args = parser.parse_args()
    # Log arguments
    logger.debug(args)
    main(args)
