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
mondo_codes = {"MONDO:0019592": "MONDO:0002145"}

# read clinvar data file from argument
df = pd.read_excel(args.input_file)
logging.info(f"Read input file: {args.input_file}")
# initial number of variants
initial_count = len(df)
logging.info(f"Initial number of variants: {initial_count}")


# filter out invalid variants

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
    f"Number of duplicate variants removed: {initial_count - len(df)}"
)
logging.info(
    f"Number of variants after filtering missing data and duplicates: {len(df)}"
)


## filter out variants that were not reported posetive
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

indel_output = (
    os.path.basename(args.input_file).replace(".xlsx","") + "_indels_50nt_or_larger.xlsx"
)
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
df["mondo_pheno"] = df["mondo_pheno"].apply(lambda x: mondo_codes.get(x, x))
logging.info("Replaced known obsolete MONDO codes")

## reofrmat date last modified column
df["LastModifiedDate"] = df["LastModifiedDate"].str.replace(
    "T00:00:00Z", "", regex=False
)
logging.info("Reformatted LastModifiedDate column")

## replace HPO terms seperator , with ;
df["Proband_HPO_terms"] = df["Proband_HPO_terms"].str.replace(
    ",", ";", regex=False
)

# if HPO:0000006 is present in Proband_HPO_terms column change "Inheritance" column to "Autosomal dominant inheritance"
mask = df["Proband_HPO_terms"].str.contains("HP:0000006")
inheritance = np.asarray(df["Inheritance"])
inheritance[mask] = "Autosomal dominant inheritance"
df["Inheritance"] = inheritance

# remove HPO:0000006 from Proband_HPO_terms column
df["Proband_HPO_terms"] = df["Proband_HPO_terms"].str.replace(
    "HP:0000006;|HP:0000006$", "", regex=True
)   

logging.info("Reformatted Proband_HPO_terms column")

# add uuid as first column
df.insert(0, "UUID", [f"uid_{uuid.uuid1().time}" for _ in range(len(df))])

# check df has no empty values otside fo the Stop and Varait_type columns
if df.drop(columns=["Stop", "Variant_type"]).isnull().values.any():
    logging.error(
        "Dataframe contains empty values after filtering and reformatting"
    )
    raise ValueError(
        "Dataframe contains empty values after filtering and reformatting"
    )
logging.info("No empty values in dataframe after filtering and reformatting")

# TODO: Remove "chr" from chromosome column if present


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
output_filtered = output_base + "_filtered.xlsx"
output_b37 = output_base + "_b37_filtered.xlsx"
output_b38 = output_base + "_b38_filtered.xlsx"
df.to_excel(output_filtered, index=False)
df_b37.to_excel(output_b37, index=False)
df_b38.to_excel(output_b38, index=False)
logging.info(f"Filtered variants written to file: {output_filtered}")
logging.info(f"GRCh37 variants written to file: {output_b37}")
logging.info(f"GRCh38 variants written to file: {output_b38}")
logging.info(f"Indels >= 50nt written to file: {indel_output}")
logging.info("Variant filtering script completed successfully")
