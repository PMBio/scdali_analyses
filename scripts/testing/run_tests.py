"""Run statistical tests for all crosses."""


import sys
import os
import argparse

import numpy as np
import scanpy as sc
import pandas as pd

from statsmodels.stats.multitest import multipletests

sys.path.append('../..')
from utils import settings

from dali import run_dali


EXP_IDS = ['SS148', 'SS157', 'SS158', 'SS159']
CELL_STATE_CONT = 'X_vae'
CELL_STATE_DISC = 'leiden_vae'
CELL_STATE_TIME = 'time_vae'
LINEAGES = ['muscle', 'nervous']
TIMEPOINTS = settings.TIMEPOINTS

REGIONS = 'cusanovich_dm6_peaks_1kb'
CORRECTION = 'wasp'
LABEL = 'f1_all_windows'


adata_path = lambda exp_id: os.path.join(
        settings.DATA_DIR,
        LABEL,
        '_'.join([exp_id, REGIONS, CORRECTION, 'allelic_counts.h5ad']))


def arg_parse():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--n_cores',
        help='Number of cores for test.',
        default=1, type=int)
    parser.add_argument(
        '--test_run',
        help='Perform test run.',
        default=False, action='store_true')
    args = parser.parse_args()
    return args


def process_adata(file_in, file_out, n_cores, test_run=False):
    adata = sc.read(file_in)

    if test_run:
        adata = adata[:, np.random.choice(adata.shape[1], 250, replace=False)].copy()
    print('[f1 analysis][TEST] Processing %d regions from %s' % (adata.shape[1], file_in))

    # count data
    A = adata.X.A.astype(float)
    D = adata.layers['allelic_total'].A.astype(float)
    autosomes_ids = (adata.var['chr'] != 'chrX').to_numpy()

    # cell state variables
    time = adata.obs[CELL_STATE_TIME].to_numpy().astype(float)
    cell_state_time = np.stack([time, time ** 2, time ** 3]).T
    cell_state_disc = pd.get_dummies(adata.obs[CELL_STATE_DISC])
    cell_state_cont = adata.obsm[CELL_STATE_CONT].astype(float)


    # run tests

    # 1. dalihet on vae embedding
    print('[f1 analysis][TEST] Running DaliHet (VAE)')
    adata.var['DALIHET_VAE'] = run_dali(
        A=A, D=D,
        cell_state=cell_state_cont,
        test='daliBB',
        n_cores=n_cores)['pvalues']

    # 2. dalihet on leiden clusters
    print('[f1 analysis][TEST] Running DaliHet (LEIDEN)')
    adata.var['DALIHET_LEIDEN'] = run_dali(
        A=A, D=D,
        cell_state=cell_state_disc,
        test='daliBB',
        n_cores=n_cores)['pvalues']

    # 3. dalihet on time
    print('[f1 analysis][TEST] Running DaliHet (TIME)')
    adata.var['DALIHET_TIME'] = run_dali(
        A=A, D=D,
        cell_state=cell_state_time,
        test='daliBB',
        n_cores=n_cores)['pvalues']

    # note: assumes no sex chromosomes present!
    # 4. dalijoint
    print('[f1 analysis][TEST] Running DaliJoint (VAE)')
    results = run_dali(
        A=A, D=D,
        cell_state=cell_state_cont,
        test='daliBB',
        base_rate=.5,
        test_cell_state_only=False,
        return_rho=True,
        n_cores=n_cores)
    adata.var['DALIJOINT'] = results['pvalues']
    adata.var['DALIJOINT_RHO'] = results['rho']

    # 5. lrt
    print('[f1 analysis][TEST] Running LRT')
    adata.var['DALIHOM'] = run_dali(
        A=A, D=D,
        test='meanBB',
        base_rate=.5,
        n_cores=n_cores)['pvalues']

    # time-lineage test
    print('[f1 analysis][TEST] Running DaliHet (TIME-LINEAGE)')
    for lineage in LINEAGES:
        lineage_cells = adata.obs['lineage_%s' % lineage].to_numpy()
        lineage_sites = adata.var['lineage_%s_covered' % lineage].to_numpy()
        # 6. lineage dalihet (time)
        pvals = run_dali(
            A=A[lineage_cells, :][:, lineage_sites],
            D=D[lineage_cells, :][:, lineage_sites],
            cell_state=cell_state_time[lineage_cells, :],
            test='daliBB',
            n_cores=n_cores)['pvalues']
        adata.var['DALIHET_TIME_%s' % lineage] = np.nan
        adata.var.loc[lineage_sites, 'DALIHET_TIME_%s' % lineage] = pvals

    if not test_run:
        adata.write(file_out)


def main():
    args = arg_parse()
    for exp_id in EXP_IDS:
        adata_file = adata_path(exp_id)
        process_adata(
                file_in=adata_file,
                file_out=adata_file,
                n_cores=args.n_cores,
                test_run=args.test_run)
        if args.test_run:
            break


if __name__ == '__main__':
    main()

