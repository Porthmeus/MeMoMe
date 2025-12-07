#!/usr/bin/env python3
"""
Helper script to recreate the tiny MYb11 submodel used in tests.
It configures COBRApy to avoid multiprocessing (required in our sandbox),
then invokes tinyModelCreator with the curated settings stored in
utils/script_accessories/tinyMyb11Settings.
"""

from __future__ import annotations

import sys
from pathlib import Path

import cobra

from utils import tinyModelCreator


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    settings_dir = project_root / "utils" / "script_accessories" / "tinyMyb11Settings"

    cobra.Configuration().processes = 1

    sys.argv = [
        tinyModelCreator.__file__,
        "--model_file_path",
        str(project_root / "tests" / "dat" / "manually_merged_models" / "MYb11_iCEL1314" / "M1_MYb11.xml"),
        "--biomass_id",
        "bio1",
        "--desired_biomass_components_file",
        str(settings_dir / "desired_biomass_components.txt"),
        "--required_excretions_file",
        str(settings_dir / "required_excretions.txt"),
        "--desired_medium_file",
        str(settings_dir / "desired_medium.csv"),
        "--keep_identifiers",
        "seed",
        "bigg",
        "chebi",
        "--new_model_id",
        "tiny_myb11",
        "--new_model_name",
        "tiny_myb11",
        "--output_file",
        str(project_root / "tests" / "dat" / "tiny_myb11.xml"),
    ]
    tinyModelCreator.main()


if __name__ == "__main__":
    main()
