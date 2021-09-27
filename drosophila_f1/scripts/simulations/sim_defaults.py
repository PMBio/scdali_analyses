import os
import sys

sys.path.append('../..')
from utils import settings


CORRECTION = 'wasp'
REGIONS = 'cusanovich_dm6_peaks_1kb'
LABEL = 'f1_all_windows'
EXP_ID = 'SS148'
ADATA_PATH = os.path.join(
    settings.DATA_DIR,
    LABEL,
    '_'.join([EXP_ID, REGIONS, CORRECTION, 'allelic_counts.h5ad']))

DEFAULT_PARAMS = {
    'test': 'DaliHet',
    'test_cellstate': 'default', # cont, disc or default (test default)
    'test_cellstate_d': 'all', # all or integer
    'ncells': 5000,
    'npeaks': 1000,
    'tcounts': 'real', # real, fixed or poisson
    'tcounts_mean': 0.1,
    'sim_kern_var': 0.1,
    'sim_overdisp': 2,
    'sim_kern_mix': 0.5,
    'sim_kern_cont': 'linear',
    'sim_kern_disc': 'global',
    'seed': 123
}

CHOICES = {
    'test': [
        'DaliJoint',
        'DaliHet',
        'DaliHom',
        'OLS',
        'BetaBinomialGLM',
        'ANOVA'],
    'test_cellstate': ['disc', 'cont', 'default'],
    'sim_kern_cont': ['RBF', 'linear'],
    'sim_kern_disc': ['global', 'cluster']
}
