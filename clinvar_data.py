import pandas as pd
import uuid
import time
import logging
import os.path
import numpy as np


class clinvar_data:

    def __init__(self, df, input_file, output_dir):
        self.df = df
        self.base_name = os.path.splitext(os.path.basename(input_file))[0]
        self.output_dir = output_dir

    @staticmethod
    def export(df, output_dir, suffix, base_name, index=False):
        output = os.path.join(output_dir, f"{base_name}{suffix}")
        df.to_excel(output, index=index)

    def filter_missing(self, df):

        #  extract invalid variants due to missing data

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

    def filter_na(self, df):
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
        """
        filter out variants that were not reported with a specific status in the summary status column.
        Args:
            df (pd.DataFrame): input dataframe
            status (str): status to filter out
        """
        df = df[df["Summary_status"] != status]
        logging.info(
            f"Number of variants after filtering out {status}: {len(df)}"
        )
        return df

    def remove_where_column_equals(self, df, column, value):
        """filter out variants where a specific column contains a specific string.
        Args:
            df (pd.DataFrame): input dataframe
            column (str): column to check for string
            string (str): string to filter out
        """
        df = df[df[column] != value]
        logging.info(
            f"Number of variants after filtering out where {column} equals '{value}': {len(df)}"
        )
        return df

    @staticmethod
    def infer_cnv_copy_number(sex, chromosome, variant_type):
        """infer copy number from the variant type, chromosome and sample sex.
        Args:
            sex (Str): Proband sex of the sample.
            chromosome (Str): Chromosome of the variant.
            variant_type (Str): Type of the variant.
        Returns:
            copy number (int) or None if cannot be inferred.
        """

        if chromosome not in ["X", "Y"]:
            if variant_type == "copy number loss":
                return "1-2"
            elif variant_type == "copy number gain":
                return "3-2"
            else:
                return None

        if sex == "FEMALE":
            if chromosome == "X":
                if variant_type == "copy number loss":
                    return "1-2"
                elif variant_type == "copy number gain":
                    return "3-2"
            if chromosome == "Y":
                return None
        if sex == "MALE":
            if chromosome == "Y":
                if variant_type == "copy number loss":
                    return "0-1"
                elif variant_type == "copy number gain":
                    return "2-1"
            if chromosome == "X":
                if variant_type == "copy number loss":
                    return "0-1"
                elif variant_type == "copy number gain":
                    return "2-1"
        return None

    def retrieve_large_variant_types(self, df, min_size: int, types: list):
        """filter varaints of particular types which are >= a minimum size.
        Args:
            df (pd.DataFrame): input dataframe
            min_size (int): minimum size of the indel
            types (list): list of variant types to filter
        """
        df_indels = df[df["Variant_type"].isin(types)]

        df_indels[["copy number", "expected copy number"]] = (
            df_indels[["Proband_sex", "Chromosome", "Variant_type"]]
            .apply(
                lambda x: self.infer_cnv_copy_number(x[0], x[1], x[2]), axis=1
            )
            .str.split("-", expand=True)
        )

        if df_indels["copy number"].isnull().any():
            logging.warning("Some copy numbers could not be inferred")

        # df_indels["copy_number"] = copy_nums

        df_indels_large = df_indels[
            (df_indels["Stop"] - df_indels["Start"]) >= min_size
        ]
        logging.info(
            f"Number of {types} larger than {min_size}nt: {len(df_indels_large)}"
        )
        return df_indels_large

    def reformat_variant_type(self, df, mapping=None, default=True):
        """Reformat the Variant_type column based on a provided mapping.

        Args:
            df (pd.DataFrame): main dataframe
            mapping (dict): dictionary mapping old variant types to new variant types
            default (bool, optional): If True, default mapping of CNV types is used. Defaults to True.
        Returns:
            pd.DataFrame: dataframe with reformatted Variant_type column
        """
        if mapping is None and default:
            mapping = {
                "deletion": "copy number loss",
                "amplification": "copy number gain",
            }

        df["Variant_type"] = df["Variant_type"].apply(
            lambda x: mapping.get(x, x)
        )

        return df

    def drop_subset(self, df, subset):
        """Remove a subset of variants from the main dataframe.

        Args:
            df (pd.DataFrame): main dataframe
            subset (pd.DataFrame): subset to remove from main dataframe

        Returns:
            pd.DataFrame: dataframe with subset removed
        """
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
        """_summary_

        Args:
            df df.DataFrame: main dataframe
            replace (bool, optional): If True, replace values in the existing column.
                If False, create a new column with reformatted values. Defaults to True.
            mapping (dict, optional): Dictionary mapping old values to new values. Defaults to None.
            column (str, optional): Name of the column to reformat. Defaults to None.
            new_column (str, optional): Name of the new column to create if replace is False. Defaults to None.
            exhaustive (bool, optional): If True, all values must be in the mapping dictionary.

        Raises:
            ValueError: If invalid parameters are provided such as both 'replace' and 'new_column' being set.

        Returns:
            pd.DataFrame: Dataframe with reformatted column.
        """
        if (
            replace
            and mapping is not None
            and column is not None
            and new_column is None
        ):
            if not exhaustive:
                df[column] = df[column].apply(
                    lambda x: mapping.get(x, "unknown")
                )
            else:
                df[column] = df[column].apply(lambda x: mapping.get(x, x))
        elif (
            new_column is not None
            and mapping is not None
            and column is not None
            and not replace
        ):
            if not exhaustive:
                df[new_column] = df[column].apply(
                    lambda x: mapping.get(x, "unknown")
                )
            else:
                df[new_column] = df[column].apply(lambda x: mapping.get(x, x))
        else:
            logging.error("Invalid parameters for reformat_columns method")
            raise ValueError("Invalid parameters for reformat_columns method")
        return df

    def replace_in_column(self, df, column, to_replace, replacement):
        """_summary_

        Args:
            df              pd.DataFrame: main dataframe
            column          str: column to perform replacement in
            to_replace      str: string to replace
            replacement     str: string to replace with

        Returns:
            pd.DataFrame : dataframe with replaced values in specified column
        """
        df[column] = df[column].str.replace(
            to_replace, replacement, regex=False
        )
        return df

    def conditional_reformat_contains(
        self, column, contains, change_col, new_value
    ):
        """change values in a column based on whether another column contains a specific string.

        Args:
            column      str: column to check for string on which conditional change is based
            contains    str: string to check for in column
            change_col  str: column to change value in in case condition is met
            new_value   str: new value to set in change_col if condition is met
        """
        mask = self.df[column].str.contains(contains, na=False)
        count = mask.sum()
        self.df.loc[mask, change_col] = new_value
        logging.info(f"{column} changed to '{new_value}' in {count} variants")

    @staticmethod
    def remove_hpo_term(df, term):
        """_remove a specific HPO term from the Proband_HPO_terms column.

        Args:
            df pd.DataFrame: input dataframe
            term str: HPO term to remove

        Returns:
            pd.DataFrame: dataframe with HPO term removed
        """
        df["Proband_HPO_terms"] = df["Proband_HPO_terms"].str.replace(
            f"{term};|{term}$", "", regex=True
        )
        logging.info(f"'{term}' removed in Proband_HPO_terms column")
        return df

    @staticmethod
    def reformat_chromosome(df):
        """reformat chromosome column by removing 'chr' prefix.

        Args:
            df df.DataFrame: input dataframe

        Returns:
            pd.DataFrame: dataframe with reformatted chromosome column
        """
        df["Chromosome"] = df["Chromosome"].apply(
            lambda x: str(x).replace("chr", "")
        )
        return df

    @staticmethod
    def insert_uuid(df):
        """insert a UUID column at the start of the dataframe.

        Args:
            df pd.DataFrame: input dataframe

        Returns:
            pd.DataFrame: dataframe with UUID column inserted
        """
        df.insert(
            0,
            "UUID",
            [
                f"uid_{uuid.uuid1().time}"
                for _ in range(len(df))
                if time.sleep(0.001) is None
            ],
        )
        return df
