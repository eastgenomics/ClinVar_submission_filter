# ClinVar_submission_filter

A script to take variant data exported from 100k Genomes and filter/clean before submission to ClinVar.

## Inputs

- ```-i```, ```--input_file``` (required): Path to an input .xlsx file of variants exported from 100k Genomes
- ```-o```, ```--output_dir``` (optional, default ```filtered```): Path to output directory to store created files in

## Example usage:
```python filter_100k_variants.py -i clinvar_data_east_of_england_prod_20250411_UNLOCKED.xlsx```

## Outputs

The script performs the following functions:
- Insert unique identifiers
- Drop rows where values for key fields are missing
- Drop rows with the Mondo generic rare disease code ```MONDO:0021136```
- Drop rows with invalid values for ```Summary_status``` or ```Classification```
- Drop duplicates
- Reformat variant classifications, CNV variant types, and obsolete Mondo codes
- Clean ```Proband_HPO_terms``` and ```LastModifiedDate``` columns
- Split variants into SNVs and CNVs and restrict by size
- Split into GRCh37 and GRCh38
- Export to .xlsx

The script produces up to 5 output .xlsx files plus a log:
- GRCh37 SNVs: ```<output_dir>/<filename>_GRCh37_SNV.xlsx```
- GRCh38 SNVs: ```<output_dir>/<filename>_GRCh38_SNV.xlsx```
- GRCh37 CNVs: ```<output_dir>/<filename>_GRCh37_CNV.xlsx```
- GRCh38 CNVs: ```<output_dir>/<filename>_GRCh38_CNV.xlsx```
- Dropped variants: ```<output_dir>/<filename>_dropped.xlsx```

# Audit script

A separate audit script has been added which can be used to compare ClinVar
reports of submitted variants against the outputs of filter_100k_variants.py.

## Requirements

- Path to outputs of filter_100k_variants.py (filtered_dir)
- Path to downloaded ClinVar submission reports (submitted_dir)

The ClinVar submission reports can be downloaded from the submission portal
and are named in the following format: ```SUBXXXXXXXX_(100)_submitter_report_B.txt```

A list of EGLH submissions can be found at https://cuhbioinformatics.atlassian.net/wiki/spaces/DV/pages/3979444228/100K+variant+submission+EGLH#%F0%9F%9A%80-List-of-submissions

## Inputs

- ```-f```, ``` --filtered ``` (required): Path to folder containing outputs of filter_100k_variants.py
- ```-s```, ``` --submitted ``` (required): Path to folder containing ClinVar submission reports

## Example usage
```python clinvar_100k_submission_audit.py -f filtered -s clinvar_reports```

## Outputs

The script prints numbers of variants which:
- Should have been submitted
- Should have been dropped
- Were actually submitted
- Were correctly submitted
- Were incorrectly submitted
- Were incorrectly not submitted
