import glob
import pandas as pd
import numpy as np

def plmc_processing(infile_loc):
    df_plmc = pd.read_csv(infile_loc, sep=' ', header=None) 
    plmc_cols = ['aa1_loc', 'trash1', 'aa2_loc', 'trash2', 'trash3', 'couplings']
    df_plmc.columns = plmc_cols
    df_plmc.drop(['trash1', 'trash2', 'trash3'], axis=1, inplace=True)
    assert all(df_plmc['aa1_loc'] < df_plmc['aa2_loc'])
    return df_plmc

def ccmpred_processing(infile_loc):
    df_ccm_pivot = pd.read_csv(infile_loc, sep='\t', header=None)
    df_ccm_pivot.drop(max(df_ccm_pivot.columns), axis=1, inplace=True)
    df_ccm_pivot.values[[np.arange(max(df_ccm_pivot.columns)+1)]*2] = np.nan
    df_ccm = df_ccm_pivot.where(np.triu(np.ones(df_ccm_pivot.shape)).astype(np.bool))
    df_ccm= df_ccm.stack().reset_index()
    df_ccm.reset_index(drop=True, inplace=True)
    df_ccm.columns = ['aa1_loc', 'aa2_loc', 'couplings']
    df_ccm['aa1_loc'] = df_ccm['aa1_loc'] + 1
    df_ccm['aa2_loc'] = df_ccm['aa2_loc'] + 1
    assert all(df_ccm['aa1_loc'] < df_ccm['aa2_loc'])
    return df_ccm

def psicov_processing(infile_loc):
    psi_cols = ['aa1_loc', 'aa2_loc', 'trash1', 'trash2', 'couplings']
    df_psi = pd.read_csv(infile_loc, sep=' ', header=None)
    df_psi.columns = psi_cols
    df_psi.sort_values(['aa1_loc','aa2_loc'], inplace=True)
    df_psi.drop(['trash1', 'trash2'], axis=1, inplace=True)
    assert all(df_psi['aa1_loc'] < df_psi['aa2_loc'])
    df_psi['aa1_loc'] = df_psi['aa1_loc'] 
    df_psi['aa2_loc'] = df_psi['aa2_loc']
    return df_psi

for infile in glob.glob('../Data/Simulated_couplings/raw/*.couplings'):
#for infile in glob.glob('../Data/Empirical_couplings/raw/*.couplings'):
    if infile.split('.')[-2] == 'ccmpred':
        df = ccmpred_processing(infile)
    
    elif infile.split('.')[-2] == 'psicov':
        try:
            temp = pd.read_csv(infile, sep=' ', header=None)    
        except pd.errors.EmptyDataError:
            continue
        df = psicov_processing(infile)
    elif infile.split('.')[-2] == 'plmc':
        df = plmc_processing(infile)
    else:
        print('ERROR')

    df.to_csv(infile.replace('.couplings', '.processed.couplings').replace('raw/', ''), index=False)

