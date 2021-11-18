import sys
import argparse

from functools import partial

import scanpy as sc
import numpy as np
import pandas as pd

import scdali

from scdali.utils.simulate import *
from scdali.utils.stats import compute_expected_sample_variance
from scdali.utils.run_model import *
from scdali.utils.parallel import process_parallel

from sim_defaults import DEFAULT_PARAMS, CHOICES, ADATA_PATH


JIT = 1e-4


#===============================================================================
# Set parameters
#===============================================================================
# check if run from snakemake
SNAKEMODE = 'snakemake' in locals() or 'snakemake' in globals()

def update_variable(var_name, default):
    """If run from snakemake, update parameter values."""
    if not SNAKEMODE:
        return default
    try:
        x = snakemake.params.simulation[var_name]
        if isinstance(x, pd.Series):
            x = x[0]
        return type(default)(x) # bit hacky
    except KeyError:
        if default is None:
            raise ValueError('%s needs to be specified' % var_name)
        else:
            return default


print('Setting parameters ...')
params = {}
for key, value in DEFAULT_PARAMS.items():
    # load default parameters & update if specified through snakemake
    params[key] = update_variable(key, value)
    try:
        if not params[key] in CHOICES[key]:
            raise ValueError('Invalid value for %s' % key)
    except KeyError:
        pass
    print('%20s: %s' % (key, str(params[key])))

params['ncores'] = snakemake.threads if SNAKEMODE else 1
params['outfile'] = snakemake.output[0] if SNAKEMODE else 'pvals.txt'

#===============================================================================
# Run simulation
#===============================================================================
print('Simulating ...')
rng = np.random.default_rng(seed=params['seed'])

# read anndata with counts and cell states
adata = sc.read(ADATA_PATH)

if params['tcounts'] == 'real':
    # filter
    adata = adata[:, adata.layers['allelic_total'].mean(0) > params['tcounts_mean']]

# subsample
if params['ncells'] > adata.shape[0] or params['npeaks'] > adata.shape[1]:
    raise ValueError('Requested simulated data size is larger than adata.')

# try not to lose any clusters, so that there are sufficient cellstate
# dimensions if cellstate is discrete
frac = params['ncells'] / adata.shape[0]
cell_ids = adata.obs.groupby('leiden_vae').sample(
    frac=frac, replace=False, random_state=params['seed']).index
adata = adata[cell_ids, :]
adata = adata[:, rng.choice(adata.shape[1], params['npeaks'], replace=False)]

if params['tcounts'] == 'real':
    D = adata.layers['allelic_total'].A
elif params['tcounts'] == 'fixed':
    D = params['tcounts_mean'] * np.ones(shape=adata.shape)
else:
    D = rng.poisson(lam=params['tcounts_mean'], size=adata.shape)
D = D.astype(int)

n, p = adata.shape

if params['sim_kern_disc'] == 'global':
    # constant kernel
    K1 = np.ones((n, n))
else:
    # cluster kernel
    cluster_ids = adata.obs['leiden_vae'].to_numpy().astype(int)
    K1 = create_cluster_kernel(cluster_ids)

# create continuous kernel
if params['sim_kern_cont'] == 'linear':
    K2 = create_linear_kernel(adata.obsm['X_vae'])
else:
    K2 = create_rbf_kernel(adata.obsm['X_vae'])

# ensure positive definiteness
jitter = JIT * np.eye(n)
K1 += jitter
K2 += jitter

# normalize kernel variance
if params['sim_kern_disc'] == 'cluster':
    K1 = K1 / compute_expected_sample_variance(K1)
K2 = K2 / compute_expected_sample_variance(K2)

# mixture kernel
K = params['sim_kern_mix'] * K1 + (1 - params['sim_kern_mix']) * K2

# simulate alternative counts
sim = simulate_beta_binomial(
    K=K, D=D,
    sigma2=params['sim_kern_var'],
    theta=params['sim_overdisp'],
    mu=0, seed=params['seed'])
A = sim['A']

# define test cellstate
if params['test_cellstate'] == 'default':
    params['test_cellstate'] = 'disc' if params['test'] == 'ANOVA' else 'cont'

if params['test_cellstate'] == 'cont':
    if params['test'] == 'ANOVA':
        raise ValueError('ANOVA requires test_cellstate disc.')
    cell_state = adata.obsm['X_vae']
else:
    cell_state = pd.get_dummies(adata.obs['leiden_vae']).to_numpy()
cell_state = cell_state.astype(float)

if params['test_cellstate_d'] == 'all':
    n_components = cell_state.shape[1]
else:
    n_components = int(float(params['test_cellstate_d']))

if n_components > cell_state.shape[1]:
    raise ValueError('Requested cell state dimensions: %d.'
        'Available: %d.' % (n_components, cell_state.shape[1]))

if params['test_cellstate'] == 'disc':
    cell_state = cell_state[:, :(n_components-1)]
    cell_state = np.hstack(
        [cell_state, cell_state.sum(1, keepdims=True) == 0])
else:
    cell_state = cell_state[:, :n_components]

# define model
if params['test'] == 'ANOVA':
    m = scdali.models.ScipyClusterTest
elif params['test'] == 'DaliHetBinomial':
    m = scdali.models.DaliHet
else:
    m = getattr(scdali.models, params['test'])

init_kwargs = {}
if params['test'] in ['DaliJoint', 'DaliHom']:
    init_kwargs['base_rate'] = .5
if params['test'] == 'DaliHetBinomial':
    init_kwargs['binomial'] = True
if params['test'] != 'DaliHom':
    init_kwargs['E'] = cell_state

callbacks = [create_method_callback('test')]
show_progress = False if params['ncores'] > 1 else True

f = partial(
        run_model,
        m,
        init_kwargs=init_kwargs,
        callbacks=callbacks,
        show_progress=show_progress)

# run
print('Testing ...')
results = process_parallel(
        f,
        mat_dict={'A':A, 'D':D},
        n_cores=params['ncores'])

pvalues = np.asarray([r[0] for r in results]).T
pd.DataFrame(pvalues).to_csv(params['outfile'], index=False, header=False)

