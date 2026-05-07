#!usr/bin/env python

import argparse
import logging
import os.path
import pandas as pd


def parse_args() -> argparse.Namespace:
    """ Parse CLI arguments.

    Returns:
        argparse.Namespace
    """

    parser = argparse.ArgumentParser(
        description="Filter and reformat an .xlsx file of Clinvar variant data"
    )

    parser.add_argument(
        "-i",
        "--input_file",
        type=str,
        help="Path to input Excel file containing ClinVar export data"
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        default="filtered",
        help="Path to output directory"
    )

    args = parser.parse_args()

    return args


def initialise_logging(output_dir):
    """ Set up logging and output log file.

    Args:
        output_dir (str): Path to output directory
    """

    # define format for log messages
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    # initialise logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # create handler for output file
    file_handler = logging.FileHandler(
        os.path.join(output_dir, "variant_filtering.log"), mode="w"
    )
    file_handler.setFormatter(formatter)

    # add handler to logger
    logger.handlers.clear()
    logger.addHandler(file_handler)


def main():

    args = parse_args()

    # ensure input file exists and is .xlsx
    if not os.path.exists(args.input_file):
        raise ValueError(f"Specified input does not exist: {args.input_file}")
    if not args.input_file.endswith('.xlsx'):
        raise ValueError("Input file must be in .xlsx format")

    # extract base name to construct output filenames
    base_name = args.input_file[:-5]

    # make output dir if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # initialise logging
    initialise_logging(args.output_dir)

    # read clinvar data into df
    logging.info(f"Reading input file: {args.input_file}")
    df = pd.read_excel(args.input_file)
    logging.info(f"File read into dataframe containing {len(df)} variants\n")

    logging.info("Variant filtering script completed successfully.")


if __name__ == "__main__":
    main()
