"""Run interpolation and variance decomposition."""


import argparse
import sys
import os

sys.path.append('../..')
from utils import settings

import scanpy as sc

from scdali import run_interpolation
from scdali.utils.stats import compute_quantile_diff, apply_fdr_bh


EXP_IDS = ['SS148', 'SS157', 'SS158', 'SS159']
CELL_STATE_CONT = 'X_vae'
CELL_STATE_DISC = 'annotation'
SUBSELECT_KEY = 'DALIHET_VAE_bh'
SUBSELECT_CUTOFF = .2


REGIONS = 'cusanovich_dm6_peaks_1kb'
CORRECTION = 'wasp'
LABEL = 'f1_all_windows'


def get_adata_path(exp_id, suffix=''):
    return os.path.join(
        settings.DATA_DIR, LABEL,
        '_'.join([exp_id, REGIONS, CORRECTION, 'allelic_counts%s.h5ad' % suffix]))


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

    # save to raw (further processing will only be applied to significantly variable sites)
    # adata.raw = adata

    # subset to significant regions
    ids = adata.var.query('%s < %f' % (SUBSELECT_KEY, SUBSELECT_CUTOFF)).index
    if test_run:
        ids = ids[:1]
    adata = adata[:, ids].copy()
    print('[f1 analysis][GP] Processing %d regions from %s' % (adata.shape[1], file_out))

    # extract counts
    A = adata.X.A.astype(float)
    D = adata.layers['allelic_total'].A.astype(float)

    cell_state_continuous = adata.obsm[CELL_STATE_CONT].astype(float)

    # run model to obtain smoothed rate estimates and effect sizes
    print('[f1 analysis][GP] Fitting interpolation model')
    results = run_interpolation(
        A=A, D=D, cell_state=cell_state_continuous,
        kernel='Linear',
        num_inducing=1100,
        max_iter=2800,
        n_cores=n_cores,
        return_prior_mean=True)

    adata.layers['gp_post_mean'] = results['posterior_mean']
    adata.layers['gp_post_var'] = results['posterior_var']
    adata.var['gp_mean'] = results['prior_mean']

    # compute effect size (5% quantile difference)
    adata.var['qdiff_5'] = compute_quantile_diff(
            adata.layers['gp_post_mean'], 0.05)
    adata.var['qdiff_10'] = compute_quantile_diff(
            adata.layers['gp_post_mean'], 0.10)

    # save
    if not test_run:
        adata.write(file_out)


def main():
    args = arg_parse()
    for exp_id in EXP_IDS:
        file_in = get_adata_path(exp_id)
        file_out = get_adata_path(exp_id, suffix='_processed')
        process_adata(
                file_in=file_in,
                file_out=file_out,
                n_cores=args.n_cores,
                test_run=args.test_run)
        if args.test_run:
            break


if __name__ == '__main__':
    main()

