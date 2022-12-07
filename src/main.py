"""
Main entry point of the program
"""
import argparse
from download_db import download


def main(args: argparse.Namespace):
    if args.download:
        print("Starting to download databases")
        download()


if __name__ == '__main__':
    # Specifies which arguments are accepted by the program
    parser = argparse.ArgumentParser(description='MeMoMe - Cool stuff.')
    # Specifying this tells the program to download all the databases
    parser.add_argument('--download', action='store_true', help='Download all required databases')

    args = parser.parse_args()
    main(args)


