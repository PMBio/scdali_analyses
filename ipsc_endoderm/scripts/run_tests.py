import pandas as pd
import numpy as np

import scanpy as sc
import scdali

# load data
adata = sc.read('../data/endoderm_ase_processed.h5ad')

def get_cell_states(N_PCS):
    E = adata.obsm['X_pca'][:, :N_PCS]
    return E / np.sqrt(np.var(E, 0))

A = adata.X.A
D = adata.layers['allelic_total'].A
X = pd.get_dummies(adata.obs['donor']).to_numpy()

for N_PCS in [1,2,3,4,5,10,15,20]:
    E = get_cell_states(N_PCS)

    adata.var[f'pval_pc{N_PCS}_dali'] = scdali.run_scdali(
        A=A, D=D, model='scDALI-Het', cell_state=E, n_cores=40)['pvalues']
    adata.var[f'pval_pc{N_PCS}_dalirs'] = scdali.run_scdali(
        A=A, D=D, X=X, model='scDALI-Het', cell_state=E, n_cores=40)['pvalues']

sc.write('../data/endoderm_ase_processed.h5ad', adata)
