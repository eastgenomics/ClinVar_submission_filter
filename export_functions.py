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


def create_snv_df(df, dropped, variant_length=50):
    """ Restrict input dataframe to SNVs under a specified length.

    Args:
        df (pd.DataFrame): Input dataframe
        dropped (list): List of dropped rows
        variant_length (int): Maximum allowable variant length. Defaults to 50.

    Returns:
        pd.DataFrame: SNV records with appropriate length
        list: Updated list of dropped rows
    """

    logging.info("Creating SNV dataframe")

    # create SNV df based on populated fields
    snv_df = df[
        df["Reference"].notnull() & \
        df["Alternate"].notnull() & \
        df["Stop"].isnull()
    ].copy()

    logging.info(f"{len(snv_df)} SNV records initially identified")

    # restrict based on variant length
    correct_length = snv_df[
        abs(snv_df['Reference'].str.len() - snv_df['Alternate'].str.len())
        < variant_length]

    drop_count = len(snv_df) - len(correct_length)

    if drop_count != 0:
        # append dropped rows to the list
        dropped_df = snv_df[~snv_df.index.isin(correct_length.index)]
        dropped_df = dropped_df.assign(drop_reason=f"SNV longer than {variant_length}nt")
        dropped += dropped_df.values.tolist()
        logging.info(f"Dropped {drop_count} SNVs longer than {variant_length}nt")

    logging.info(f"SNV dataframe created with {len(correct_length)} records\n")

    return correct_length, dropped


def create_cnv_df(df, dropped, variant_length=50):
    """ Restrict input dataframe to CNVs of at least a specified length.

    Args:
        df (pd.DataFrame): Input dataframe
        dropped (list): List of dropped rows
        variant_length (int): Minimum allowable variant length. Defaults to 50.

    Returns:
        pd.DataFrame: CNV records with appropriate length
        list: Updated list of dropped rows
    """

    logging.info("Creating CNV dataframe")

    # select records with CNV variant type
    cnv_df = df[df["Variant_type"].isin(
        ['copy number gain', 'copy number loss']
    )].copy()

    initial_cnvs = len(cnv_df)
    logging.info(f"{initial_cnvs} CNV records initially identified")

    # check correct fields are populated
    correct_fields = cnv_df[
        cnv_df["Reference"].isnull() & \
        cnv_df["Alternate"].isnull() & \
        cnv_df["Stop"].notnull()
    ]

    # add records with wrong fields to dropped list
    wrong_fields = len(cnv_df) - len(correct_fields)
    if wrong_fields != 0:
        dropped_df = cnv_df[~cnv_df.index.isin(correct_fields.index)]
        dropped_df = dropped_df.assign(drop_reason="CNV with incorrect fields populated")
        dropped += dropped_df.values.tolist()
        logging.info(f"Dropped {wrong_fields} CNVs with incorrect fields")

    # restrict based on variant length
    correct_length = correct_fields[(correct_fields["Stop"] - correct_fields["Start"]) >= variant_length]

    # add records below required length to dropped list
    too_short = len(correct_fields) - len(correct_length)
    if too_short != 0:
        dropped_df = correct_fields[~correct_fields.index.isin(correct_length.index)]
        dropped_df = dropped_df.assign(drop_reason=f"CNV shorter than {variant_length}nt")
        dropped += dropped_df.values.tolist()
        logging.info(f"Dropped {too_short} CNVs shorter than {variant_length}nt")

    # infer CNV copy numbers
    correct_length[["copy number", "expected copy number"]] = correct_length.apply(
        lambda x: infer_cnv_copy_number(
            x.Proband_sex, x.Chromosome, x.Variant_type),
            axis=1, result_type='expand')

    if correct_length["copy number"].isnull().any():
        logging.warning("Some copy numbers could not be inferred")

    logging.info(f"CNV dataframe created with {len(correct_length)} records\n")

    return correct_length, dropped


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
