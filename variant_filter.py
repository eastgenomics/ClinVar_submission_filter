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
    obsolete_codes := "--obsolete-codes", type=str, help="Path to the obsolete MONDO codes mapping file or string of paired values", default=None
)

## parse the arguments
args = parser.parse_args()

# setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    filename=args.output_dir + "/variant_filter.log",
    filemode="w",
)
logging.basicConfig(
    level=logging.ERROR,
    format="%(asctime)s - %(levelname)s - %(message)s",
    filename=args.output_dir + "/variant_filter.err",
    filemode="w",
)
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

# obselete MONDO codes mapping
## check if obsolete codes file is provided or string of paired values
mondo_codes = None
if args.obsolete_codes:
    if os.path.isfile(args.obsolete_codes):
        try:
            mondo_df = pd.read_csv(args.obsolete_codes, sep="\t", header=None)
            mondo_codes = dict(zip(mondo_df[0], mondo_df[1]))
            logging.info(f"Read obsolete MONDO codes from file: {args.obsolete_codes}")
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
            logging.info("Read obsolete MONDO codes from string of paired values")
        except Exception as e:
            logging.error(f"Error parsing obsolete MONDO codes string: {e}")
            raise


# read clinvar data file from argument
df = pd.read_excel(args.input_file)
logging.info(f"Read input file: {args.input_file}")
# initial number of variants
initial_count = len(df)
logging.info(f"Initial number of variants: {initial_count}")


# filter out invalid variants

## write variants with missing data to separate file
df_missing = df[df[["Start",
        "Chromosome",
        "mondo_pheno",
        "Classification",
        "Build"]].isnull().any(axis=1)]

## filter out variants with missing data
df = df.dropna(
    subset=[
        "Start",
        "Chromosome",
        "mondo_pheno",
        "Classification",
        "Build",
    ]
)
logging.info(
    f"Number of variants with missing data removed: {initial_count - len(df)}"
)
rolling_count = len(df)

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
logging.info(
    f"Number of duplicate variants removed: {rolling_count - len(df)}"
)
rolling_count = len(df)

logging.info(
    f"Number of variants after filtering missing data and duplicates: {len(df)}"
)


## filter out variants that were not reported  positive
df = df[df["Summary_status"] != "REPORTED_INCONCLUSIVE"]
logging.info(f"Number of variants after filtering out REPORTED_INCONCLUSIVE: {len(df)}")
logging.info(f"variants filtered out: {initial_count - len(df)}")

## filter indels >= 50nt to separate file
df_indels = df[
    (df["Variant_type"] == "deletion")
    | (df["Variant_type"] == "amplification")
]
logging.info(f"Number of indels: {len(df_indels)}")

df_indels_large = df_indels[(df_indels["Stop"] - df_indels["Start"]) >= 50]

base_name = os.path.splitext(os.path.basename(args.input_file))[0]
indel_output = os.path.join(args.output_dir, f"{base_name}_indels_50nt_or_larger.xlsx")

df_indels_large.to_excel(indel_output, index=False)

logging.info(f"Number of indels >= 50nt: {len(df_indels_large)}")
logging.info(f"Indels >= 50nt written to file: {indel_output}")

## drop indels >= 50nt from main dataframe
df = df[~df.index.isin(df_indels_large.index)]
logging.info(f"Number of variants after removing indels >= 50nt: {len(df)}")

# reformat
## reformat clinical significance column
df["Classification_reformated"] = df["Classification"].apply(
    lambda x: clinsigs.get(x, "Unknown")
)
logging.info("Reformatted clinical significance column")
logging.info(
    f"Number of variants with unknown clinical significance: \
        {len(df[df['Classification_reformated'] == 'Unknown'])}"
)


## replace known invaid MONDO codes
if mondo_codes:
    df["Mondo_code"] = df["Mondo_code"].apply(lambda x: mondo_codes.get(x, x))
    logging.info("Replaced known obsolete MONDO codes")

## reformat date last modified column
df["LastModifiedDate"] = df["LastModifiedDate"].str.replace(
    "T00:00:00Z", "", regex=False
)
logging.info("Reformatted LastModifiedDate column")

## replace HPO terms separator , with ;
df["Proband_HPO_terms"] = df["Proband_HPO_terms"].str.replace(
    ",", ";", regex=False
)

# if HPO:0000006 is present in Proband_HPO_terms column change "Inheritance" column to "Autosomal dominant inheritance"

mask = df["Proband_HPO_terms"].str.contains("HP:0000006", na=False)
df.loc[mask, "Inheritance"] = "Autosomal dominant inheritance"

# remove HPO:0000006 from Proband_HPO_terms column
df["Proband_HPO_terms"] = df["Proband_HPO_terms"].str.replace(
    "HP:0000006;|HP:0000006$", "", regex=True
)   

logging.info("Reformatted Proband_HPO_terms column")


# Remove "chr" from chromosome column if present
df["Chromosome"] = df["Chromosome"].apply(lambda x: str(x).replace("chr", ""))

# add uuid as first column
df.insert(0, "UUID", [f"uid_{uuid.uuid1().time}" for _ in range(len(df))])

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
output_base = args.output_dir + os.path.basename(args.input_file).replace(".xlsx", "")
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
