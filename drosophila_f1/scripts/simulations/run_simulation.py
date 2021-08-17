import sys
import argparse

import scanpy as sc
import numpy as np
import pandas as pd

from dali.utils.simulate import *
from dali.utils.stats import compute_expected_sample_variance
from dali import run_dali


JIT = 1e-4


def parse_args():
    parser = argparse.ArgumentParser()

    # model to run, model kernel, data, outfile
    parser.add_argument('--model',
            choices=['dalihet', 'lrt', 'dalijoint', 'ttest'],
            help='Test to run.',
            required=True)
    parser.add_argument('--adata',
            help='Anndata file with counts and cell state.',
            required=True)
    parser.add_argument('--outprefix',
            help='Output file for p-values.',
            required=True)

    # simulation parameters
    parser.add_argument('--sim_kern_var',
            help='Kernel variance component.',
            type=float, required=True)
    parser.add_argument('--sim_overdisp',
            help='Overdispersion parameter.',
            type=float, required=True)
    parser.add_argument('--sim_kern_mix',
            help=('Mixture coefficient to interpolate between continuous'
                'and cluster kernel. If 0, use continous, if 1, use discrete.'),
            type=float, required=True)
    parser.add_argument('--sim_kern_cont',
            choices=['RBF', 'linear'],
            help='Type of continuous kernel.',
            default='linear')
    parser.add_argument('--sim_kern_disc',
            choices=['cluster', 'global'],
            help='Type of continuous kernel.')


    # misc
    parser.add_argument('--ncores',
            help='Number of cores.',
            type=int, default=1)
    parser.add_argument('--seed',
            help='Random seed.',
            type=int, default=123)
    parser.add_argument('--test',
            help='Run test.',
            default=False,
            action='store_true')

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    # read anndata with counts and cell states
    adata = sc.read(args.adata)
    if args.test:
        adata = adata[:, np.random.choice(adata.shape[1], 50)]
    n, p = adata.shape

    print('define kernel')
    if args.sim_kern_disc == 'global':
        # constant kernel
        K1 = np.ones((n, n))
    else:
        # cluster kernel
        cluster_ids = adata.obs['leiden_vae'].to_numpy().astype(int)
        K1 = create_cluster_kernel(cluster_ids)

    # create continuous kernel
    if args.sim_kern_cont == 'linear':
        K2 = create_linear_kernel(adata.obsm['X_vae'])
    else:
        K2 = create_rbf_kernel(adata.obsm['X_vae'])

    # ensure positive definiteness
    jitter = JIT * np.eye(n)
    K1 += jitter
    K2 += jitter

    # normalize kernel variance TODO: only if not global!
    if args.sim_kern_disc == 'cluster':
        K1 = K1 / compute_expected_sample_variance(K1)
    K2 = K2 / compute_expected_sample_variance(K2)

    # mixture kernel
    K = args.sim_kern_mix * K1 + (1 - args.sim_kern_mix) * K2
    D = adata.layers['allelic_total'].A.astype(int)

    print('simulate')
    # simulate alternative counts
    sim = simulate_beta_binomial(
        K=K, D=D,
        sigma2=args.sim_kern_var,
        theta=args.sim_overdisp,
        mu=0, seed=args.seed)
    A = sim['A']

    if args.model == 'dalihet':
        results = run_dali(
            A=A, D=D,
            cell_state=adata.obsm['X_vae'].astype(float),
            test='daliBB',
            n_cores=args.ncores)
    elif args.model == 'dalijoint':
        results = run_dali(
            A=A, D=D,
            cell_state=adata.obsm['X_vae'].astype(float),
            test='daliBB',
            test_cell_state_only=False,
            return_rho=True,
            base_rate=.5,
            n_cores=args.ncores)
    elif args.model == 'lrt':
        results = run_dali(
            A=A, D=D,
            test='meanBB',
            base_rate=.5,
            n_cores=args.ncores)
    elif args.model == 'ttest':
        results = run_dali(
            A=A, D=D,
            cell_state=pd.get_dummies(adata.obs['leiden_vae']).to_numpy().astype(float),
            test='ttest',
            n_cores=args.ncores)

    if args.model == 'dalijoint':
        pd.DataFrame(results['rho']).to_csv(args.outprefix + '-result_rhos.txt')
    pd.DataFrame(results['pvalues']).to_csv(args.outprefix + '-result_pvals.txt')


if __name__ == '__main__':
    main()
