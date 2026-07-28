# ClinVar_submission_filter

A script to take variant data exported from 100k Genomes and filter/clean before submission to ClinVar.

## Inputs

- ```-i```, ```--input_file``` (required): Path to an input .xlsx file of variants exported from 100k Genomes
- ```-o```, ```--output_dir``` (optional, default ```filtered```): Path to output directory to store created files in

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

The script produces 5 output .xlsx files plus a log:
- GRCh37 SNVs: ```<output_dir>/<filename>_GRCh37_SNV.xlsx```
- GRCh38 SNVs: ```<output_dir>/<filename>_GRCh38_SNV.xlsx```
- GRCh37 CNVs: ```<output_dir>/<filename>_GRCh37_CNV.xlsx```
- GRCh38 CNVs: ```<output_dir>/<filename>_GRCh38_CNV.xlsx```
- Dropped variants: ```<output_dir>/<filename>_dropped.xlsx```
