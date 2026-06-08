#!/usr/bin/env python

import logging
import pandas as pd


def drop_missing_data_rows(df, dropped, required_fields):
    """ Drop rows where required fields are missing.

    Args:
        df (pd.DataFrame): Input dataframe
        dropped (list): List of dropped rows
        required_fields (list): Fields which must be populated, specified as
            REQUIRED_FIELDS in config.py

    Returns:
        pd.DataFrame: Dataframe after dropping rows with missing values
        list: Updated list of dropped rows
    """

    logging.info("Dropping rows where key fields are missing")

    before_drop = len(df)

    # drop rows where required fields are NA
    retained_df = df.dropna(subset=required_fields, how="any")

    # append dropped rows to the list
    dropped_df = df[~df.index.isin(retained_df.index)]
    dropped_df = dropped_df.assign(drop_reason="Missing required field(s)")
    dropped += dropped_df.values.tolist()

    logging.info(f"Number of records dropped: {before_drop - len(retained_df)}\n")

    return retained_df, dropped


def drop_duplicates(df, dropped, key_fields):
    """ Drop duplicates based on a specified combination of fields.

    Args:
        df (pd.DataFrame): Input dataframe
        dropped (list): List of dropped rows
        key_fields (list): Fields used to determine duplicate rows, specified
            as DUPLICATE_FIELDS in config.py

    Returns:
        pd.DataFrame: Dataframe after dropping duplicates
        list: Updated list of dropped rows
    """

    logging.info("Dropping duplicate records")

    before_drop = len(df)

    # ensure column to sort by is in datetime format
    df["LastModifiedDate"] = pd.to_datetime(df["LastModifiedDate"], errors='coerce')

    # drop duplicate rows
    retained_df = df.sort_values("LastModifiedDate", na_position='first').drop_duplicates(subset=key_fields, keep='last')

    # append dropped rows to the list
    dropped_df = df[~df.index.isin(retained_df.index)]
    dropped_df = dropped_df.assign(drop_reason="Duplicate record")
    dropped += dropped_df.values.tolist()

    logging.info(f"Number of records dropped: {before_drop - len(retained_df)}\n")

    return retained_df, dropped


def drop_generic_mondo(df, dropped, generic_code):
    """ Drop rows where the Mondo code is the generic code for rare disease.

    Args:
        df (pd.DataFrame): Input dataframe
        dropped (list): List of dropped rows
        generic_code (str): Mondo code for generic rare disease indication,
            specified in config.py as GENERIC_MONDO

    Returns:
        pd.DataFrame: Dataframe after removing rows
        list: Updated list of dropped rows
    """

    logging.info("Dropping records with generic Mondo code for rare disease")

    before_drop = len(df)

    # drop records with generic mondo codes
    retained_df = df[df["Mondo_code"] != generic_code]

    # append dropped rows to the list
    dropped_df = df[~df.index.isin(retained_df.index)]
    dropped_df = dropped_df.assign(drop_reason="Mondo code for generic RD indication")
    dropped += dropped_df.values.tolist()

    logging.info(f"Number of records dropped: {before_drop - len(retained_df)}\n")

    return retained_df, dropped


def drop_where_column_equals(df, dropped, column, value):
    """ Remove rows where specified column has specified value.

    Args:
        df (pd.DataFrame): Input dataframe
        dropped (list): List of dropped rows
        column (str): Column to check for string
        value (str): Value in column to filter out variants.

    Returns:
        pd.DataFrame: Dataframe after removing rows
        list: Updated list of dropped rows
    """

    logging.info(f"Dropping records where column {column} has value {value}")

    before_drop = len(df)

    # remove rows where specified column has specified value
    retained_df = df[df[column] != value]

    # append dropped rows to the list
    dropped_df = df[~df.index.isin(retained_df.index)]
    dropped_df = dropped_df.assign(drop_reason=f"Column {column} has value {value}")
    dropped += dropped_df.values.tolist()

    logging.info(f"Number of records dropped: {before_drop - len(retained_df)}\n")

    return retained_df, dropped
