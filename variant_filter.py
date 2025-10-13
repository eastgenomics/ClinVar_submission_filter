# improt
import pandas as pd
import uuid
import time
import argparse
import logging
import os.path
import numpy as np
from clinvar_data import clinvar_data

# from clinvar_data import clinvar_data

# set arguments
## create the parser
parser = argparse.ArgumentParser(
    description="Filter and reformat Clinvar data"
)
## add the arguments
parser.add_argument(
    "input_file", type=str, help="Path to the input Clinvar Excel file"
)

parser.add_argument(
    "output_dir", type=str, help="Path to the output directory"
)

parser.add_argument(
    obsolete_codes := "--obsolete-codes",
    type=str,
    help="Path to the obsolete MONDO codes mapping file or string of paired values",
    default=None,
)

## parse the arguments
args = parser.parse_args()

# ensure output dir exists
os.makedirs(args.output_dir, exist_ok=True)

# setup logging with separate files
logger = logging.getLogger()
logger.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

info_handler = logging.FileHandler(
    os.path.join(args.output_dir, "variant_filter.log"), mode="w"
)
info_handler.setLevel(logging.INFO)
info_handler.setFormatter(formatter)

error_handler = logging.FileHandler(
    os.path.join(args.output_dir, "variant_filter.err"), mode="w"
)
error_handler.setLevel(logging.ERROR)
error_handler.setFormatter(formatter)

logger.handlers.clear()
logger.addHandler(info_handler)
logger.addHandler(error_handler)

logging.info("Starting variant filtering script")


# clinical significance mapping
clinsigs = {
    "benign_variant": "Benign",
    "likely_benign_variant": "Likely benign",
    "variant_of_unknown_clinical_significance": "Uncertain significance",
    "likely_pathogenic_variant": "Likely pathogenic",
    "pathogenic_variant": "Pathogenic",
    "not_assessed": "not provided",
}


mondo_codes = None
# obselete MONDO codes mapping
## check if obsolete codes file is provided or string of paired values
if args.obsolete_codes:
    if os.path.isfile(args.obsolete_codes):
        try:
            mondo_df = pd.read_csv(args.obsolete_codes, sep="\t", header=None)
            mondo_codes = dict(zip(mondo_df[0], mondo_df[1]))
            logging.info(
                f"Read obsolete MONDO codes from file: {args.obsolete_codes}"
            )
        except Exception as e:
            logging.error(f"Error reading obsolete MONDO codes file: {e}")
            raise
    else:
        try:
            pairs = args.obsolete_codes.split(",")
            mondo_codes = {}
            for pair in pairs:
                logging.info(f"Parsing pair: {pair}")
                old, new = pair.split("-")
                mondo_codes[old] = new
            logging.info(
                "Read obsolete MONDO codes from string of paired values"
            )
        except Exception as e:
            logging.error(f"Error parsing obsolete MONDO codes string: {e}")
            raise


# read clinvar data file from argument
df = pd.read_excel(args.input_file)
logging.info(f"Read input file: {args.input_file}")
# initial number of variants
initial_count = len(df)
logging.info(f"Initial number of variants: {initial_count}")


def main(
    df=df,
    count=initial_count,
    input_file=args.input_file,
    output_dir=args.output_dir,
):
    data = clinvar_data(df, input_file, output_dir)

    # write variants with missing data to separate file
    missing_data = data.filter_missing(data.df)
    # filtered variants with missing data
    rolling_count = len(data.df)
    logging.info(
        f"Number of variants with missing data removed: {count - len(missing_data)}"
    )
    # filter out varaints with missing data in required columns
    data.df = data.filter_na(data.df)
    # filter out duplicat variants
    data.df = data.filter_duplicates(data.df)

    logging.info(
        f"Number of duplicate variants removed: {rolling_count - len(data.df)}"
    )
    rolling_count = len(df)

    logging.info(
        f"Number of variants after filtering missing data and duplicates: {len(data.df)}"
    )

    data.df = data.remove_reported_with(
        df=data.df, status="REPORTED_INCONCLUSIVE"
    )
    logging.info(f"variants filtered out: {rolling_count - len(data.df)}")
    rolling_count = len(data.df)

    cnv_data = data.retrieve_variant_types(
        df=data.df,
        type=["amplification", "deletion", "insertion"],
        min_size=50,
    )
    logging.info(f"CNVs >= 50nt written to file in {output_dir}")

    # export filtered indels >= 50nt
    data.export(
        cnv_data,
        output_dir=args.output_dir,
        sufix="_indels_50nt.xlsx",
        index=False,
    )

    data.df = data.drop_subset(data.df, cnv_data)
    logging.info(
        f"Number of variants after removing indels >= 50nt: {len(data.df)}"
    )
    rolling_count = len(data.df)

    data.df = data.reformat_columns(
        data.df,
        replace=False,
        mapping=clinsigs,
        column="Classification",
        new_column="Classification_reformated",
        exhaustive=False,
    )
    logging.info("Reformatted clinical significance column")
    logging.info(
        f"Number of variants with unknown clinical significance: \
        {len(data.df[data.df['Classification_reformated'] == 'Unknown'])}"
    )

    data.df = data.reformat_columns(
        data.df,
        replace=True,
        mapping=clinsigs,
        column="Mondo_code",
        exhaustive=False,
    )
    logging.info("Replaced known obsolete MONDO codes")

    data.df = data.replace_in_column(
        data.df, "LastModifiedDate", "T00:00:00Z", ""
    )

    logging.info("Reformatted LastModifiedDate column")

    data.df = data.replace_in_column(data.df, "Proband_HPO_terms", ",", ";")

    logging.info("Replaced seperator ',' with ';' in Proband_HPO_terms column")

    # Add class and function calls

    # check df has no empty values outside of the Stop and Variant_type columns
    if data.df.drop(columns=["Stop", "Variant_type"]).isnull().values.any():
        logging.error(
            "Dataframe contains empty values after filtering and reformatting"
        )
        raise ValueError(
            "Dataframe contains empty values after filtering and reformatting"
        )
    logging.info(
        "No empty values in dataframe after filtering and reformatting"
    )

    # final number of variants
    final_count = len(df)
    logging.info(f"Final number of variants: {final_count}")

    # split filtered table into b37 and b38 using the "Build" column
    df_b37 = df[df["Build"] == "GRCh37"]
    df_b38 = df[df["Build"] == "GRCh38"]
    logging.info(f"Number of GRCh37 variants: {len(df_b37)}")
    logging.info(f"Number of GRCh38 variants: {len(df_b38)}")
    # write to output file
    clinvar_data.export(
        missing_data,
        output_dir,
        "_missing_data.xlsx",
        index=False,
        base_name=data.base_name,
    )
    clinvar_data.export(
        df, output_dir, "_filtered.xlsx", index=False, base_name=data.base_name
    )
    clinvar_data.export(
        df_b37,
        output_dir,
        "_b37_filtered.xlsx",
        index=False,
        base_name=data.base_name,
    )
    clinvar_data.export(
        df_b38,
        output_dir,
        "_b38_filtered.xlsx",
        index=False,
        base_name=data.base_name,
    )
    logging.info(
        f"Filtered variants written to file: {data.base_name}_filtered.xlsx"
    )
    logging.info(
        f"GRCh37 variants written to file: {data.base_name}_b37_filtered.xlsx"
    )
    logging.info(
        f"GRCh38 variants written to file: {data.base_name}_b38_filtered.xlsx"
    )
    logging.info(
        f"Indels >= 50nt written to file: {data.base_name}_indels_50nt.xlsx"
    )
    logging.info("Variant filtering script completed successfully")


if __name__ == "__main__":
    main()
