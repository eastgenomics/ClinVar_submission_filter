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

    def retrieve_large_variant_types(self, df, min_size: int, types: list):
        """filter varaints of particular types which are >= a minimum size.
        Args:
            df (pd.DataFrame): input dataframe
            min_size (int): minimum size of the indel
            types (list): list of variant types to filter
        """
        df_indels = df[df["Variant_type"].isin(types)]

        df_indels_large = df_indels[
            (df_indels["Stop"] - df_indels["Start"]) >= min_size
        ]
        logging.info(
            f"Number of {types} larger than {min_size}nt: {len(df_indels_large)}"
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
        return df

    def replace_in_column(self, df, column, to_replace, replacement):
        df[column] = df[column].str.replace(
            to_replace, replacement, regex=False
        )
        return df

    def conditional_reformat_contains(
        self, column, contains, change_col, new_value
    ):
        mask = self.df[column].str.contains(contains, na=False)
        count = mask.sum()
        self.df.loc[mask, change_col] = new_value
        logging.info(f"{column} changed to '{new_value}' in {count} variants")

    @staticmethod
    def remove_hpo_term(df, term):
        df["Proband_HPO_terms"] = df["Proband_HPO_terms"].str.replace(
            f"{term};|{term}$", "", regex=True
        )
        logging.info(f"'{term}' removed in Proband_HPO_terms column")
        return df

    @staticmethod
    def reformat_chromosome(df):
        df["Chromosome"] = df["Chromosome"].apply(
            lambda x: str(x).replace("chr", "")
        )
        return df

    @staticmethod
    def insert_uuid(df):
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
