#!/usr/bin/env python

import argparse
import logging
import os.path
import pandas as pd
import config

from row_drop_functions import (
    drop_missing_data_rows,
    drop_duplicates,
    drop_generic_mondo,
    drop_where_column_equals
)

from reformat_functions import (
    reformat_column,
    replace_single_column_value,
    insert_uuid
)

from export_functions import (
   create_snv_df,
   create_cnv_df,
   split_builds,
   export_df_to_xl
)


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
        required=True,
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
    base_name = os.path.splitext(os.path.basename(args.input_file))[0]

    # make output dir if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # initialise logging
    initialise_logging(args.output_dir)

    # read clinvar data into df
    logging.info(f"Reading input file: {args.input_file}")
    df = pd.read_excel(args.input_file)
    logging.info(f"File read into dataframe containing {len(df)} variants\n")

    # initialise list of dropped rows
    df['drop_reason'] = ""
    dropped = []

    # drop rows with missing, duplicate, or inappropriate values
    df, dropped = drop_missing_data_rows(df, dropped, config.REQUIRED_FIELDS)
    df, dropped = drop_duplicates(df, dropped, config.DUPLICATE_FIELDS)
    df, dropped = drop_generic_mondo(df, dropped, config.GENERIC_MONDO)
    df, dropped = drop_where_column_equals(df, dropped, "Summary_status", "REPORTED_INCONCLUSIVE")
    df, dropped = drop_where_column_equals(df, dropped, "Classification", "not_assessed")

    # clean data values
    df = reformat_column(df, config.MONDO_MAP, "Mondo_code")
    df = reformat_column(df, config.CLINSIG_MAP, "Classification", new_column="Classification_reformatted", replace=False, stringent=True)
    df = reformat_column(df, config.CNV_MAP, "Variant_type")
    df = replace_single_column_value(df, "LastModifiedDate", " 00:00:00+00:00", "")
    df = replace_single_column_value(df, "Proband_HPO_terms", ",", ";")

    # add UUID column
    df = insert_uuid(df)

    # filter by variant type and genome build
    snv_df, dropped = create_snv_df(df, dropped)
    cnv_df, dropped = create_cnv_df(df, dropped)

    logging.info("Splitting SNV records into GRCh37/GRCh38 dataframes")
    snv_b37, snv_b38 = split_builds(snv_df)
    logging.info("Splitting CNV records into GRCh37/GRCh38 dataframes")
    cnv_b37, cnv_b38 = split_builds(cnv_df)

    # export df of dropped variants
    drop_df = pd.DataFrame(dropped, columns=df.columns)
    drop_df = replace_single_column_value(drop_df, "LastModifiedDate", " 00:00:00+00:00", "")
    export_df_to_xl(drop_df, args.output_dir, base_name, "_dropped.xlsx")

    # export files
    export_df_to_xl(snv_b37, args.output_dir, base_name, "_GRCh37_SNV.xlsx")
    export_df_to_xl(snv_b38, args.output_dir, base_name, "_GRCh38_SNV.xlsx")
    export_df_to_xl(cnv_b37, args.output_dir, base_name, "_GRCh37_CNV.xlsx")
    export_df_to_xl(cnv_b38, args.output_dir, base_name, "_GRCh38_CNV.xlsx")

    logging.info("Variant filtering script completed successfully.")


if __name__ == "__main__":
    main()
