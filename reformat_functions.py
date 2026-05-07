#!/usr/bin/env python

import logging
import pandas as pd
import uuid
import time


def reformat_column(
        df, mapping, column, new_column=None, replace=True, stringent=False
    ) -> pd.DataFrame:

    """ Given a specified df column and a mapping of which values to
    replace, either replace column values in-place, or create a new column
    to hold replacement values.

    Args:
        df (pd.DataFrame): Input dataframe
        mapping (dict): Mapping of old values to new values
        column (str): Name of column to reformat
        new_column (str): Name of new column to create if replace=False
        replace (bool): If True, modify column in-place; if False, create new
            column to hold reformatted values
        stringent (bool): If True, all column values must be in mapping dict

    Raises:
        BOTH_SET_ERROR: If both 'replace' and 'new_column' are not None
        NEITHER_SET_ERROR: If both 'replace' and 'new_column' are None

    Returns:
        pd.DataFrame: Dataframe with new or reformatted column
    """

    logging.info(f"Reformatting column: {column}\n")

    # raise error if invalid parameters provided
    BOTH_SET_ERROR = "'replace' and 'new_column' are both set but are mutually exclusive"
    NEITHER_SET_ERROR = "'replace' and 'new_column' are both None, one must be provided"

    if replace and new_column:
        logging.error(BOTH_SET_ERROR)
        raise ValueError(BOTH_SET_ERROR)

    elif not new_column and not replace:
        logging.error(NEITHER_SET_ERROR)
        raise ValueError(NEITHER_SET_ERROR)

    if replace:
        target_column = column
    elif new_column:
        target_column = new_column

    if not stringent:
        df[target_column] = df[column].apply(
            lambda x: mapping.get(x, "unknown")
        )

        unknown = len(df[df[target_column] == 'unknown'])
        if unknown != 0:
            logging.info(f"{unknown} rows had values not present in the map")

    else:
        df[target_column] = df[column].apply(lambda x: mapping.get(x,x))

    return df


def replace_single_column_value(
        df, column, to_replace, replacement
    ) -> pd.DataFrame:

    """ Replace a specified value within a specified column.

    Args:
        df (pd.DataFrame): Input dataframe
        column (str): Column to perform replacement in
        to_replace (str): String to replace
        replacement (str): String to replace with

    Returns:
        pd.DataFrame : Dataframe with replaced values in specified column
    """

    logging.info(f"Replacing {to_replace} with {replacement} in {column} column\n")

    df[column] = df[column].astype("string")
    df[column] = df[column].str.replace(to_replace, replacement, regex=False)

    return df


def insert_uuid(df) -> pd.DataFrame:
    """ Insert a UUID column at the start of a dataframe.

    Args:
        df (pd.DataFrame): Input dataframe

    Returns:
        pd.DataFrame: Dataframe with UUID column inserted
    """

    logging.info(f"Adding UUID column to dataframe\n")

    uuids = []

    for _ in range(len(df)):
        uuids.append(f"uid_{uuid.uuid1().time}")
        time.sleep(0.001)

    df.insert(0, "UUID", uuids)

    return df
