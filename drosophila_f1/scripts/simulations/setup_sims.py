"""Create pandas dataframes with simulation parameters."""


from sim_defaults import DEFAULT_PARAMS
from itertools import product
import numpy as np
import pandas as pd

SIM_OVERDISP = [2.0, 5.0]
SIM_KERN_VAR = [0.1, 0.05, 0.01]
SIM_KERN_MIX = np.linspace(0, 1, 6)
TEST_CELLSTATE_D = np.linspace(2, 24, 12)
NCELLS = [250, 500, 1000, 5000]
TCOUNTS_MEAN = 5.0
SEED = [DEFAULT_PARAMS['seed']]
SEED = SEED + np.random.default_rng(123).integers(0, 10000, 9).tolist()

COLUMNS = DEFAULT_PARAMS.keys()

def add_default_params(df):
    for c in COLUMNS:
        if c not in df.columns:
            df[c] = DEFAULT_PARAMS[c]
    return df

df_all = None
# # ==============================================================================
# # POWER experiments
# # ==============================================================================

# # heterogeneous vs. homogeneous effects
# tests = ['DaliJoint', 'DaliHet', 'DaliHom']
# df = pd.DataFrame(product(
#     tests, SIM_OVERDISP, SIM_KERN_VAR, SIM_KERN_MIX, SEED))
# df.columns = ['test', 'sim_overdisp', 'sim_kern_var', 'sim_kern_mix', 'seed']
# df['sim_kern_disc'] = 'global'
# df.to_csv('params/params_hethom.csv', index=False)
# df_all = pd.concat([df_all, add_default_params(df)])

# # continuous vs. discrete effects
# tests = ['DaliHet', 'ANOVA', 'OLS', 'BetaBinomialGLM']
# df = pd.DataFrame(product(
#     tests, SIM_OVERDISP, SIM_KERN_VAR, SIM_KERN_MIX, SEED))
# df.columns = ['test', 'sim_overdisp', 'sim_kern_var', 'sim_kern_mix', 'seed']
# df['sim_kern_disc'] = 'cluster'
# df.to_csv('params/params_contdisc.csv', index=False)
# df_all = pd.concat([df_all, add_default_params(df)])

# # ==============================================================================
# # CALIBRATION experiments
# # ==============================================================================

# # no heterogeneous imbalance
# tests = ['DaliHet', 'ANOVA', 'OLS', 'BetaBinomialGLM']
# df = pd.DataFrame(product(
#     tests, SIM_OVERDISP, SIM_KERN_VAR))
# df.columns = ['test', 'sim_overdisp', 'sim_kern_var']
# df['sim_kern_disc'] = 'global'
# df['sim_kern_mix'] = 1.0
# df.to_csv('params/params_null1.csv', index=False)
# df_all = pd.concat([df_all, add_default_params(df)])

# # neither heterogeneous nor homogeneous imbalance
# tests = ['DaliJoint', 'DaliHom']
# df = pd.DataFrame(product(
#     tests, SIM_OVERDISP))
# df.columns = ['test', 'sim_overdisp']
# df['sim_kern_disc'] = 'global'
# df['sim_kern_var'] = 0.0
# df.to_csv('params/params_null2.csv', index=False)
# df_all = pd.concat([df_all, add_default_params(df)])

# calibration for varying numbers of ncells and test_cellstate_d
tests = ['DaliHet', 'ANOVA', 'OLS', 'BetaBinomialGLM']
df = pd.DataFrame(product(
    tests, TEST_CELLSTATE_D, NCELLS, SIM_OVERDISP, SEED))
df.columns = ['test', 'test_cellstate_d', 'ncells', 'sim_overdisp', 'seed']
df['test_cellstate'] = 'disc'
df['sim_kern_var'] = 0.0
df['sim_kern_disc'] = 'global'
df['tcounts'] = 'fixed'
df['tcounts_mean'] = TCOUNTS_MEAN
df.to_csv('params/params_null3.csv', index=False)
df_all = pd.concat([df_all, add_default_params(df)])

df_all.drop_duplicates().to_csv('params/params_all.csv', index=False)
