#!/usr/bin/env python

"""
This script performs the following actions:
- Reads in output files from ClinVar_submission_filter (i.e. which variants
should have been submitted or dropped)
- Reads in ClinVar reports detailing variants whose submission was accepted
- Identifies variants which were submitted, but should have been dropped
- Identifies variants which should have been submitted, but weren't
- Confirms that all submitted variants are from the lists of those which
should have been either submitted or dropped (i.e. we know their provenance)

Outputs:
- Counts and CSV files of variants in each category

Requirements:
- Path to outputs of ClinVar_submission_filter (filtered)
- Path to downloaded ClinVar submission reports (submitted)

The ClinVar submission reports are named in the following format:
- SUBXXXXXXXX_(100)_submitter_report_B.txt
"""

import argparse
import glob
import pandas as pd
import subprocess
import re
import os


def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-f", "--filtered", type=str, required=True,
        help="Path to folder containing outputs of filter_100k_variants.py"
    )

    parser.add_argument(
        "-s", "--submitted", type=str, required=True,
        help="Path to folder containing ClinVar submission reports"
    )

    args = parser.parse_args()

    return args


def read_vars_to_submit(dir, cols):
    """ Read the 3 ClinVar_submission_filter outputs which list the variants to
    submit. Combine their records into a single df. """

    # identify, read and concatenate the files
    to_submit_fps = glob.glob(f"{dir}/*GRCh*.xlsx")
    assert len(to_submit_fps) == 3, f"Expected 3 files, found {len(to_submit_fps)}"

    to_submit = pd.concat(
        [pd.read_excel(x, usecols=cols) for x in to_submit_fps],
        ignore_index=True
    )

    # Stop value is required for matching ClinVar records
    to_submit['Stop'] = to_submit.apply(
        lambda x: get_stop_value(x['Start'], x['Stop'], x['Reference']),
        axis=1
    )

    # get rid of ref/alt nan values
    to_submit['Reference'] = to_submit['Reference'].apply(lambda x: str(x).replace('nan', ''))
    to_submit['Alternate'] = to_submit['Alternate'].apply(lambda x: str(x).replace('nan', ''))

    # standardise column dtypes
    to_submit = set_field_types(to_submit)

    # print number of variants
    to_submit.to_csv('audit_to_submit.csv', index=False)
    print(f"Variants to submit: {len(to_submit)}")

    return to_submit


def read_vars_to_drop(dir, cols):
    """ Read the ClinVar_submission_filter output which lists variants that
    should be dropped. """

    # read in variants that should have been dropped
    to_drop_fps = glob.glob(f"{dir}/*dropped.xlsx")
    assert len(to_drop_fps) == 1, f"Expected 1 file, found {len(to_drop_fps)}"
    to_drop = pd.read_excel(to_drop_fps[0], usecols=cols)

    # can ignore variants missing these fields because ClinVar wouldn't accept them
    to_drop = to_drop.dropna(how='any',
        subset=["Build", "Chromosome", "Start", "Mondo_code", "Classification"]
    )

    # reconstruct Classification_reformatted column
    mapping = {
        'benign_variant': 'Benign',
        'likely_benign_variant': 'Likely benign',
        'variant_of_unknown_clinical_significance': 'Uncertain significance',
        'likely_pathogenic_variant': 'Likely pathogenic',
        'pathogenic_variant': 'Pathogenic',
        'not_assessed': 'not provided'
    }

    to_drop['Classification_reformatted'] = to_drop['Classification'].apply(lambda x: mapping.get(x, "unknown"))

    # Stop value is required for matching ClinVar records
    to_drop['Stop'] = to_drop.apply(
        lambda x: get_stop_value(x['Start'], x['Stop'], x['Reference']),
        axis=1
    )

    # get rid of ref/alt nan values
    to_drop['Reference'] = to_drop['Reference'].apply(lambda x: str(x).replace('nan', ''))
    to_drop['Alternate'] = to_drop['Alternate'].apply(lambda x: str(x).replace('nan', ''))

    # standardise column dtypes
    to_drop = set_field_types(to_drop)

    # print number of variants
    to_drop.to_csv('audit_to_drop.csv', index=False)
    print(f"Variants to drop: {len(to_drop)}")

    return to_drop


def read_submitted_files(dir):
    """ Read the ClinVar submission reports which list variants whose
    submissions were accepted. Combine into a single df. """

    # get paths to submission records
    submitted_fps = glob.glob(f"{dir}/*_submitter_report_B.txt")
    submitted_list = []

    # read in all submission reports
    for fp in submitted_fps:

        # remove commas so they don't mess up the output CSV
        new_fp = fp.replace('B.txt', 'B_fixed.txt')
        with open(new_fp, 'w') as writer:
            subprocess.run(["sed", "-e", "s/,//g", fp], stdout=writer)

        single_df = pd.read_csv(new_fp, sep='\t', header=25)

        # add column for which submission a record was in
        single_df['submission'] = re.findall(r'SUB\d{8}(?=_)', os.path.basename(fp))[0]
        submitted_list.append(single_df)

    submitted = pd.concat(submitted_list, ignore_index=True)

    # rename columns for readability/consistency
    submitted = submitted.rename(columns={
        'Your_record_id': 'UUID_submitted',
        'Your_variant_description_chromosome_coordinates': 'variant_info',
        'Your_condition_identifier': 'Mondo_code',
        'ClinVar_condition_name': 'ClinVar_pheno',
        'Germline_classification': 'Classification_reformatted',
        'Date_last_evaluated': 'LastModifiedDate'
    })

    # split variant info into separate columns
    submitted[[
        'Build', 'Chromosome', 'Start', 'Stop', 'Reference', 'Alternate'
    ]] = submitted['variant_info'].apply(clean_variant_string).tolist()

    # remove extra MONDO: prefix
    submitted['Mondo_code'] = submitted['Mondo_code'].apply(lambda x: str(x)[6:])

    # there are a lot of unnecessary fields
    submitted = submitted.drop(columns=[
        '#Your_variant_id', 'VariationID', 'AlleleID', 'RCV', 'variant_info',
        'Your_variant_description_HGVS', 'Preferred_variant_name',
        'Your_condition_name', 'Assigned_Concept_ID', 'Assertion_criteria',
        'Submission_date', 'Submitted_gene', 'Novelty', 'Status', 'Info'
    ])

    # standardise column dtypes
    submitted = set_field_types(submitted)

    # print number of variants
    submitted.to_csv('audit_submitted.csv', index=False)
    print(f"Variants submitted: {len(submitted)}")

    return submitted


def get_stop_value(start, stop, ref):
    """ If Stop is available, convert it to integer. If not, identify it based
    on Start and length of the Reference allele. """
    try:
        return int(stop)
    except ValueError:
        return int(start + (len(str(ref)) - 1))


def clean_variant_string(var_string):
    """ Split variant details from submission data into separate columns.

    Unfortunately some records don't contain Reference or Alternate values, so
    we can't just use split() because list positions aren't consistent... """

    # get rid of quotes
    var_string = var_string.replace('"', '')

    # get fields which are reliably present
    build = re.findall(r'(?<=Assembly=)(.*?)(?=\s)', var_string)[0]
    chrom = re.findall(r'(?<=Chr=)(.*?)(?=\s)', var_string)[0]
    start = re.findall(r'(?<=start=)(.*?)(?=\s)', var_string)[0]
    stop = re.findall(r'(?<=stop=)(.*?)(?=\s)', var_string)[0]

    # ref and alt may or may not be there
    try:
        ref = re.findall(r'(?<=referenceAllele=)(.*?)(?=\s)', var_string)[0]
    except IndexError:
        ref = ""
    try:
        alt = re.findall(r'(?<=alternateAllele=)(.*?)(?=\s)', var_string)[0]
    except IndexError:
        alt = ""

    output = [build, chrom, start, stop, ref, alt]

    return output


def set_field_types(df):
    """ Standardise column dtypes to enable merging later on """

    # Start & Stop fields should be integers
    df = df.astype({'Start': int, 'Stop': int}, errors='ignore')

    # date field should be YYYY-MM-DD
    df['LastModifiedDate'] = df['LastModifiedDate'].apply(lambda x: x.replace('T00:00:00Z', ''))
    df['LastModifiedDate'] = pd.to_datetime(df['LastModifiedDate'], yearfirst=True)

    # all other columns should be strings
    for col in df.columns:
        if col not in ['Start', 'Stop', 'LastModifiedDate']:
            df[col] = df[col].apply(lambda x: str(x).strip())

    return df


def get_incorrectly_not_submitted(to_submit, submitted):
    """ Identify records from the list of variants to submit which were either
    correctly submitted, or incorrectly not submitted. Match on as many fields
    as possible. """

    # get df of records correctly submitted
    correct = to_submit.merge(submitted, how='inner', on=[
        'Build', 'Chromosome', 'Start', 'Stop', 'Reference', 'Alternate',
        'Mondo_code', 'Classification_reformatted', 'LastModifiedDate'
    ])

    correct = correct[[
        'Build', 'Chromosome', 'Start', 'Stop', 'Reference', 'Alternate',
        'Mondo_code', 'Classification_reformatted', 'LastModifiedDate',
        'Summary_status', 'SCV', 'UUID_submitted', 'submission'
    ]]

    # get df of records incorrectly not submitted (left vs inner merge)
    missing = to_submit.merge(submitted, how='left', indicator=True, on=[
        'Build', 'Chromosome', 'Start', 'Stop', 'Reference', 'Alternate',
        'Mondo_code', 'Classification_reformatted', 'LastModifiedDate'
    ])

    # drop records which appear in the df of submitted variants
    missing = missing[missing['_merge'] == 'left_only'].drop(columns=['_merge'])

    missing = missing[[
            'submission', 'SCV', 'UUID_submitted', 'Summary_status', 'Build',
            'Chromosome', 'Start', 'Stop', 'Reference', 'Alternate', 'Mondo_code',
            'Classification', 'Classification_reformatted', 'LastModifiedDate'
        ]]

    correct.to_csv('audit_correctly_submitted.csv', index=False)
    missing.to_csv('audit_incorrectly_not_submitted.csv', index=False)

    print(f"Variants correctly submitted: {len(correct)}")
    print(f"Variants incorrectly not submitted: {len(missing)}")

    return correct, missing


def get_incorrectly_submitted(to_drop, submitted, correctly_submitted):
    """ Get records which should have been dropped, but appear in the df of
    submitted variants. Match on as many fields as possible. """

    # inner merge keeps only records present in both dfs
    df = to_drop.merge(submitted, how='inner', on=[
        'Build', 'Chromosome', 'Start', 'Stop', 'Reference', 'Alternate',
        'Mondo_code', 'Classification_reformatted', 'LastModifiedDate'
    ])

    # restrict to relevant columns
    df = df[[
        'Build', 'Chromosome', 'Start', 'Stop', 'Reference', 'Alternate',
        'Mondo_code', 'Classification_reformatted', 'LastModifiedDate',
        'Summary_status', 'SCV', 'UUID_submitted', 'submission', 'drop_reason'
    ]]

    # some dropped variants are exact duplicates of correctly submitted ones
    # so it is correct for those variants to have been submitted
    df = df[~df['SCV'].isin(correctly_submitted['SCV'])]

    df.to_csv('audit_incorrectly_submitted.csv', index=False)
    print(f"Variants incorrectly submitted: {len(df)}")

    return df


def get_mystery_variants(to_submit, to_drop, submitted):
    """ Confirm that all submitted variants appeared in either the to_submit
    or to_drop dataframes (i.e. we know their provenance). Return any which
    didn't. """

    cols = [
        'Build', 'Chromosome', 'Start', 'Stop', 'Reference', 'Alternate',
        'Mondo_code', 'Classification_reformatted', 'LastModifiedDate'
    ]

    # exclude records present in to_submit
    df1 = submitted.merge(to_submit, how='left', indicator=True, on=cols)
    df1 = df1[df1['_merge'] == 'left_only'].drop(columns=['_merge'])

    # exclude records present in to_drop
    df2 = df1.merge(to_drop, how='left', indicator=True, on=cols)
    df2 = df2[df2['_merge'] == 'left_only'].drop(columns=['_merge'])

    # print any remaining
    df2.to_csv('audit_mystery_variants.csv', index=False)
    print(f"Variants submitted with unclear provenance: {len(df2)}")

    return df2


def main():

    args = parse_args()

    # read ClinVar_submission_filter outputs
    filtered_cols = [
        'Build', 'Chromosome', 'Start', 'Stop', 'Reference', 'Alternate',
        'Summary_status', 'Mondo_code', 'LastModifiedDate', 'Classification',
        'Classification_reformatted', 'drop_reason'
    ]

    to_submit = read_vars_to_submit(args.filtered, filtered_cols)
    to_drop = read_vars_to_drop(args.filtered, filtered_cols)

    # read ClinVar submission reports (accepted variants)
    submitted = read_submitted_files(args.submitted)

    # identify variants which should have been submitted but weren't
    correctly_submitted, incorrectly_not_submitted = get_incorrectly_not_submitted(to_submit, submitted)

    # identify variants which were submitted but shouldn't have been
    incorrectly_submitted = get_incorrectly_submitted(to_drop, submitted, correctly_submitted)

    # confirm all submitted variants were from the filtering outputs
    left_over = get_mystery_variants(to_submit, to_drop, submitted)

    # the numbers in stdout don't add up because some incorrectly submitted
    # variants were submitted MULTIPLE TIMES
    check = incorrectly_submitted.drop_duplicates(
        subset=[
            'Build', 'Chromosome', 'Start', 'Stop', 'Reference', 'Alternate',
            'Mondo_code', 'Classification_reformatted', 'LastModifiedDate'
        ]
    )

    # confirm this results in the numbers adding up correctly
    assert len(submitted) == (len(correctly_submitted) + len(check) + len(left_over))
    assert len(to_submit) == (len(correctly_submitted) + len(incorrectly_not_submitted))


if __name__ == '__main__':
    main()
