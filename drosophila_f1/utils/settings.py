"""
================================================================================
Analysis settings, constants and directory structures.
================================================================================
"""

import os
#===============================================================================
# Files used in the analysis. Needs to be modified!
#===============================================================================
# peak files
TOTAL_COUNTS_PEAKS = '/icgc/dkfzlsdf/analysis/B260/users/heinent/projects/f1_pipeline/data/resources/cusanovich_dm6_peaks.txt'
ALLELIC_COUNTS_PEAKS = '/icgc/dkfzlsdf/analysis/B260/users/heinent/projects/f1_pipeline/data/resources/cusanovich_dm6_peaks_1kb.txt'
CUSANOVICH_PEAKS = '/icgc/dkfzlsdf/analysis/B260/users/heinent/projects/f1_pipeline/data/resources/cusanovich_dm6_peaks_unsorted.txt'

# cell barcodes files
CELL_INDEX = lambda exp_id: f'/icgc/dkfzlsdf/analysis/B260/users/heinent/projects/f1_pipeline/data/processed/{exp_id}/cell_barcodes/{exp_id}.readdepth_auto.cells_indextable.txt'

# counts files
TOTAL_COUNTS = lambda exp_id: f'/icgc/dkfzlsdf/analysis/B260/users/heinent/projects/f1_pipeline/data/processed/{exp_id}/counts/{exp_id}.cusanovich_dm6_peaks.nowasp.counts_total.txt.gz'
ALLELE1_COUNTS = lambda exp_id: f'/icgc/dkfzlsdf/analysis/B260/users/heinent/projects/f1_pipeline/data/processed/{exp_id}/counts/{exp_id}.cusanovich_dm6_peaks_1kb.wasp.counts_allele1.txt.gz'
ALLELE2_COUNTS = lambda exp_id: f'/icgc/dkfzlsdf/analysis/B260/users/heinent/projects/f1_pipeline/data/processed/{exp_id}/counts/{exp_id}.cusanovich_dm6_peaks_1kb.wasp.counts_allele2.txt.gz'

# variant info files
VARIANT_INFO_LENIENT = lambda exp_id: f'/icgc/dkfzlsdf/analysis/B260/users/heinent/projects/f1_pipeline/data/processed/{exp_id}/vcf_filtered/{exp_id}.variant_info_lenient.txt.gz'
VARIANT_INFO_STRINGENT = lambda exp_id: f'/icgc/dkfzlsdf/analysis/B260/users/heinent/projects/f1_pipeline/data/processed/{exp_id}/vcf_filtered/{exp_id}.variant_info_stringent.txt.gz'
#===============================================================================
# Constants
#===============================================================================
# timepoint-stage correspondence
TIMEPOINT_TO_STAGE = {
    '2to4' : ['stage4-6', 'stage7-8'],
    '6to8' : ['stage9-10', 'stage11-12'],
    '10to12' : ['stage13-16'],
}

# available timepoints
TIMEPOINTS = [
    '2to4',
    '6to8',
    '10to12'
]

# experiment ids and parental strains
EXP_IDS_TO_PATERNAL_DGRP = {
    'SS148': 'DGRP-712',
    'SS157': 'DGRP-639',
    'SS158': 'DGRP-852',
    'SS159': 'DGRP-307'
}
F1_EXP_IDS = EXP_IDS_TO_PATERNAL_DGRP.keys()


#===============================================================================
# Directory structure
#===============================================================================
ANALYSIS_ROOT_DIR = '/'.join(os.path.abspath(__file__).split('/')[:-2])
DATA_DIR = os.path.join(ANALYSIS_ROOT_DIR, 'data')
RESOURCES_DIR = os.path.join(DATA_DIR, 'resources')

