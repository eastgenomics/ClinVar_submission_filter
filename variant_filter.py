# improt
import pandas as pd
import uuid
import time
import argparse
import logging
import os.path
import numpy as np

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
    missing_data = data.filter_missing()
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

    data = data.drop_subset(data.df, cnv_data)
    logging.info(
        f"Number of variants after removing indels >= 50nt: {len(data.df)}"
    )
    rolling_count = len(data.df)

    data.df = data.reformat_columns(
        data.df,
        replace=True,
        mapping=clinsigs,
        column="Classification",
        new_column="Classification_reformated",
        exhaustive=False,
    )
    logging.info("Reformatted clinical significance column")
    logging.info(
        f"Number of variants with unknown clinical significance: \
        {len(df[df['Classification_reformated'] == 'Unknown'])}"
    )

    data.df = data.reformat_columns(
        data.df,
        replace=True,
        mapping=clinsigs,
        column="mondo_code",
        new_column=mondo_codes,
        exhaustive=False,
    )
    logging.info("Replaced known obsolete MONDO codes")

    data = data.replace_in_column(
        data.df, "LastModifiedDate", "T00:00:00Z", ""
    )

    logging.info("Reformatted LastModifiedDate column")

    data = data.replace_in_column(data.df, "Proband_HPO_terms", ",", ";")

    logging.info("Replaced seperator ',' with ';' in Proband_HPO_terms column")


# TODO: finish refactoring below code into methods of clinvar_data class, them move to its own module
class clinvar_data:

    def __init__(
        self, df, input_file=args.input_file, output_dir=args.output_dir
    ):
        self.df = df
        self.base_name = os.path.splitext(os.path.basename(input_file))[0]

    @staticmethod
    def export(self, df, output_dir, sufix, index=False):
        output = os.path.join(output_dir, f"{self.base_name}{sufix}")
        df.to_excel(output, index=index)

    def filter_missing(self, df):

        # filter out invalid variants

        ## write variants with missing data to separate file
        filter_df = df[
            df[
                [
                    "Start",
                    "Chromosome",
                    "mondo_pheno",
                    "Classification",
                    "Build",
                ]
            ]
            .isnull()
            .any(axis=1)
        ]
        return filter_df

    def filter_na(df):
        ## filter out variants with missing data
        filter_df = df.dropna(
            subset=[
                "Start",
                "Chromosome",
                "mondo_pheno",
                "Classification",
                "Build",
            ]
        )
        return filter_df

    def filter_duplicates(self, df):
        ## filter out duplicates
        df = df.drop_duplicates(
            subset=[
                "Start",
                "Stop",
                "Chromosome",
                "Reference",
                "Alternate",
                "mondo_pheno",
            ],
            keep="first",
        )
        return df

    def remove_reported_with(self, df, status):
        ## filter out variants that were not reported  positive
        df = df[df["Summary_status"] != status]
        logging.info(
            f"Number of variants after filtering out {status}: {len(df)}"
        )
        return df

    def retrieve_variant_types(
        self, df, min_size: int = None, type: list = None
    ):
        ## filter indels >= 50nt to separate file
        df_indels = df[(df["Variant_type"] in type)]

        df_indels_large = df_indels[
            (df_indels["Stop"] - df_indels["Start"]) >= 50
        ]
        logging.info(
            f"Number of {type} larger than {min_size}nt: {len(df_indels_large)}"
        )
        return df_indels_large

    def drop_subset(self, df, subset):
        df = df[~df.index.isin(subset.index)]
        return df

    def reformat_columns(
        self,
        df,
        replace=True,
        mapping=None,
        column=None,
        new_column=None,
        exhaustive=False,
    ):
        if replace and mapping and column and not new_column:
            if not exhaustive:
                df[column] = df[column].apply(
                    lambda x: mapping.get(x, "unknown")
                )
            else:
                df[column] = df[column].apply(lambda x: mapping.get(x, x))
        elif new_column and mapping and column and not replace:
            if not exhaustive:
                df[new_column] = df[column].apply(
                    lambda x: mapping.get(x, "unknown")
                )
            else:
                df[new_column] = df[column].apply(lambda x: mapping.get(x, x))
        else:
            logging.error("Invalid parameters for reformat_columns method")
            raise ValueError("Invalid parameters for reformat_columns method")

    def replace_in_column(self, df, column, to_replace, replacement):
        df[column] = df[column].str.replace(
            to_replace, replacement, regex=False
        )
        return df


# TODO: refactor below code into functions and methods of clinvar_data class

# if HPO:0000006 is present in Proband_HPO_terms column change "Inheritance" column to "Autosomal dominant inheritance"

mask = df["Proband_HPO_terms"].str.contains("HP:0000006", na=False)
adi_count = mask.sum()
df.loc[mask, "Inheritance"] = "Autosomal dominant inheritance"

logging.info(
    f"Inheritance changed to 'Autosomal dominant inheritance' in {adi_count} variants"
)

# remove HPO:0000006 from Proband_HPO_terms column
df["Proband_HPO_terms"] = df["Proband_HPO_terms"].str.replace(
    "HP:0000006;|HP:0000006$", "", regex=True
)

logging.info("'HP:0000006' removed in Proband_HPO_terms column")


# Remove "chr" from chromosome column if present
df["Chromosome"] = df["Chromosome"].apply(lambda x: str(x).replace("chr", ""))

# add uuid as first column
df.insert(
    0,
    "UUID",
    [
        f"uid_{uuid.uuid1().time}"
        for _ in range(len(df))
        if time.sleep(0.001) is None
    ],
)

# check df has no empty values outside of the Stop and Variant_type columns
if df.drop(columns=["Stop", "Variant_type"]).isnull().values.any():
    logging.error(
        "Dataframe contains empty values after filtering and reformatting"
    )
    raise ValueError(
        "Dataframe contains empty values after filtering and reformatting"
    )
logging.info("No empty values in dataframe after filtering and reformatting")


# final number of variants
final_count = len(df)
logging.info(f"Final number of variants: {final_count}")

# split filtered table into b37 and b38 using the "Build" column
df_b37 = df[df["Build"] == "GRCh37"]
df_b38 = df[df["Build"] == "GRCh38"]
logging.info(f"Number of GRCh37 variants: {len(df_b37)}")
logging.info(f"Number of GRCh38 variants: {len(df_b38)}")
# write to output file
base_name = os.path.splitext(os.path.basename(args.input_file))[0]
output_base = os.path.join(args.output_dir, base_name)
output_missing = output_base + "_missing_data.xlsx"
output_filtered = output_base + "_filtered.xlsx"
output_b37 = output_base + "_b37_filtered.xlsx"
output_b38 = output_base + "_b38_filtered.xlsx"
df_missing.to_excel(output_missing, index=False)
df.to_excel(output_filtered, index=False)
df_b37.to_excel(output_b37, index=False)
df_b38.to_excel(output_b38, index=False)
logging.info(f"Filtered variants written to file: {output_filtered}")
logging.info(f"GRCh37 variants written to file: {output_b37}")
logging.info(f"GRCh38 variants written to file: {output_b38}")
logging.info(f"Indels >= 50nt written to file: {indel_output}")
logging.info("Variant filtering script completed successfully")


if __name__ == "__main__":
    main()
