#!/usr/bin/env python

import logging
import pandas as pd


def drop_missing_data_rows(df, required_fields) -> pd.DataFrame:
    """ Drop rows where required fields are missing.

    Args:
        df (pd.DataFrame): Input dataframe
        required_fields (list): Fields which must be populated, specified as
            REQUIRED_FIELDS in config.py

    Returns:
        pd.DataFrame: Dataframe after dropping rows with missing values
    """

    logging.info("Dropping rows where key fields are missing")

    before_drop = len(df)

    df.dropna(subset=required_fields, how="any", inplace=True)

    logging.info(f"Number of records dropped: {before_drop - len(df)}\n")

    return df


def drop_duplicates(df, key_fields) -> pd.DataFrame:
    """ Drop duplicates based on a specified combination of fields.

    Args:
        df (pd.DataFrame): Input dataframe
        key_fields (list): Fields used to determine duplicate rows, specified
            as DUPLICATE_FIELDS in config.py

    Returns:
        pd.DataFrame: Dataframe after dropping duplicates
    """

    logging.info("Dropping duplicate records")

    before_drop = len(df)

    # ensure column to sort by is in datetime format
    df["LastModifiedDate"] = pd.to_datetime(df["LastModifiedDate"], errors='coerce')

    df = df.sort_values("LastModifiedDate", na_position='first').drop_duplicates(subset=key_fields, keep='last')

    logging.info(f"Number of records dropped: {before_drop - len(df)}\n")

    return df


def drop_generic_mondo(df, generic_code) -> pd.DataFrame:
    """ Drop rows where the Mondo code is the generic code for rare disease.

    Args:
        df (pd.DataFrame): Input dataframe
        generic_code (str): Mondo code for generic rare disease indication,
            specified in config.py as GENERIC_MONDO

    Returns:
        pd.DataFrame: Dataframe after removing rows
    """

    logging.info("Dropping records with generic Mondo code for rare disease")

    before_drop = len(df)

    df = df[df["Mondo_code"] != generic_code]

    logging.info(f"Number of records dropped: {before_drop - len(df)}\n")

    return df


def drop_where_column_equals(df, column, value) -> pd.DataFrame:
    """ Remove rows where specified column has specified value.

    Args:
        df (pd.DataFrame): Input dataframe
        column (str): Column to check for string
        value (str): Value in column to filter out variants.

    Returns:
        pd.DataFrame: Dataframe after removing rows
    """

    logging.info(f"Dropping records where column {column} has value {value}")

    before_drop = len(df)

    df = df[df[column] != value]

    logging.info(f"Number of records dropped: {before_drop - len(df)}\n")

    return df
