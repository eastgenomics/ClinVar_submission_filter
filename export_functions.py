#!/usr/bin/env python

import logging
import os.path
import pandas as pd


def infer_cnv_copy_number(sex, chromosome, variant_type) -> list:
    """ Infer copy number (in a very general sense) based on variant type,
    chromosome and sample sex.

    Args:
        sex (str): Proband sex of the sample
        chromosome (str): Chromosome of the variant
        variant_type (str): Type of variant (copy number gain or loss)

    Returns:
        list or None: "copy number-expected copy number" where expected copy
            number is based on chromosome and sex. Returns None if values
            cannot be inferred.
    """

    if sex == 'FEMALE' and chromosome == 'Y':
        return None

    # handle autosomes, or X in females
    if chromosome not in ["X", "Y"] \
        or ((sex == "FEMALE") and (chromosome == 'X')):

        if variant_type == 'copy number gain':
            return ['3', '2']

        elif variant_type == 'copy number loss':
            return ['1', '2']

    # handle sex chromosome variants in males
    elif chromosome in ["X", "Y"] and sex == 'MALE':

        if variant_type == 'copy number gain':
            return ['2', '1']

        elif variant_type == 'copy number loss':
            return ['0', '1']

    return None


def create_snv_df(df, variant_length=50) -> pd.DataFrame:
    """ Restrict input dataframe to SNVs under a specified length.

    Args:
        df (pd.DataFrame): Input dataframe
        variant_length (int): Maximum allowable variant length. Defaults to 50.

    Returns:
        pd.DataFrame: SNV records with appropriate length
    """

    logging.info("Creating SNV dataframe")

    # create SNV df based on populated fields
    snv_df = df[
        df["Reference"].notnull() & \
        df["Alternate"].notnull() & \
        df["Stop"].isnull()
    ].copy()

    initial_snvs = len(snv_df)
    logging.info(f"SNV dataframe created with {initial_snvs} records")

    # restrict based on variant length
    snv_df = snv_df[
        abs(snv_df['Reference'].str.len() - snv_df['Alternate'].str.len())
        < variant_length]

    dropped = initial_snvs - len(snv_df)
    if dropped != 0:
        logging.info(f"Dropped {dropped} SNVs longer than {variant_length}nt\n")

    return snv_df


def create_cnv_df(df, variant_length=50) -> pd.DataFrame:
    """ Restrict input dataframe to CNVs of at least a specified length.

    Args:
        df (pd.DataFrame): Input dataframe
        variant_length (int): Minimum allowable variant length. Defaults to 50.

    Returns:
        pd.DataFrame: CNV records with appropriate length
    """

    logging.info("Creating CNV dataframe")

    # create CNV df based on populated fields
    cnv_df = df[
        df["Reference"].isnull() & \
        df["Alternate"].isnull() & \
        df["Stop"].notnull()
    ].copy()

    initial_cnvs = len(cnv_df)
    logging.info(f"CNV dataframe created with {initial_cnvs} records\n")

    # restrict based on variant length
    cnv_df = cnv_df[(cnv_df["Stop"] - cnv_df["Start"]) >= variant_length]

    long_cnvs = len(cnv_df)
    too_short = initial_cnvs - long_cnvs
    if too_short != 0:
        logging.info(f"Dropped {too_short} CNVs shorter than {variant_length}nt")

    # ensure CNV records have an appropriate variant type
    cnv_df = cnv_df[cnv_df["Variant_type"].isin(
        ['copy number gain', 'copy number loss']
    )]

    wrong_type = long_cnvs - len(cnv_df)
    if wrong_type != 0:
        logging.info(f"Dropped {wrong_type} CNVs with wrong variant type")

    # infer CNV copy numbers
    cnv_df[["copy number", "expected copy number"]] = cnv_df.apply(
        lambda x: infer_cnv_copy_number(
            x.Proband_sex, x.Chromosome, x.Variant_type),
            axis=1, result_type='expand')

    if cnv_df["copy number"].isnull().any():
        logging.warning("Some copy numbers could not be inferred")

    return cnv_df


def split_builds(df) -> pd.DataFrame:
    """ Given a variant dataframe, split GRCh37/GRCh38 variants into separate
    dataframe.

    Args:
        df (pd.DataFrame): Input dataframe

    Returns:
        pd.DataFrame: Two separate GRCh37/GRCh38 dataframes
    """

    b37_df = df[df["Build"] == "GRCh37"]
    b38_df = df[df["Build"] == "GRCh38"]

    if len(b37_df) + len(b38_df) != len(df):
        BUILD_ERROR = "Some build values are not GRCh37/GRCh38"
        logging.error(BUILD_ERROR)
        raise ValueError(BUILD_ERROR)

    # sort dfs
    fields = [
        "Chromosome", "Start", "Stop", "Reference", "Alternate", "Classification"
    ]

    b37_df = b37_df.sort_values(fields, ignore_index=True)
    b38_df = b38_df.sort_values(fields, ignore_index=True)

    # log final variant numbers
    logging.info(f"Number of GRCh37 variants: {len(b37_df)}")
    logging.info(f"Number of GRCh38 variants: {len(b38_df)}\n")

    return b37_df, b38_df


def export_df_to_xl(df, output_dir, base_name, suffix):
    """ Export a Pandas df to an Excel file.

    Args:
        df (pd.dataframe): Dataframe to export
        output_dir (str): Path to store exported file under
        base_name (str): Prefix for exported file
        suffix (str): Suffix for exported file
    """

    output_path = os.path.join(output_dir, f"{base_name}{suffix}")

    if len(df) == 0:
        logging.warning(f"Empty dataframe, no export to {output_path}")

    else:
        logging.info(f"Exporting {len(df)} records to {output_path}")
        df.to_excel(output_path, index=False)
