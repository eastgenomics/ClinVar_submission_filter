REQUIRED_FIELDS = [
    "Build",
    "Chromosome",
    "Start",
    "Mondo_code",
    "Classification"
]

DUPLICATE_FIELDS = [
    "Build",
    "Chromosome",
    "Start",
    "Stop",
    "Reference",
    "Alternate",
    "Variant_type",
    "mondo_pheno",
    "Classification"
]

GENERIC_MONDO = "MONDO:0021136"

MONDO_MAP = {
    'MONDO:0019592': 'MONDO:0002145',
    'MONDO:0018734': 'MONDO:0009076'
}

CLINSIG_MAP = {
    'benign_variant': 'Benign',
    'likely_benign_variant': 'Likely benign',
    'variant_of_unknown_clinical_significance': 'Uncertain significance',
    'likely_pathogenic_variant': 'Likely pathogenic',
    'pathogenic_variant': 'Pathogenic',
    'not_assessed': 'not provided'
}

CNV_MAP = {
    'deletion': 'copy number loss',
    'amplification': 'copy number gain'
}
