import glob
import pandas as pd
import numpy as np

def processing(infile_loc):
    psi_cols = ['aa1_loc', 'aa2_loc', 'trash1', 'trash2', 'couplings']
    df_psi = pd.read_csv(infile_loc, sep=' ', header=None)
    df_psi.columns = psi_cols
    df_psi.sort_values(['aa1_loc','aa2_loc'], inplace=True)
    df_psi.drop(['trash1', 'trash2'], axis=1, inplace=True)
    assert all(df_psi['aa1_loc'] < df_psi['aa2_loc'])
    df_psi['aa1_loc'] = df_psi['aa1_loc'] 
    df_psi['aa2_loc'] = df_psi['aa2_loc']
    return df_psi

for infile in glob.glob('../Data/Empirical_ml/raw/*.couplings'):
    try:
        temp = pd.read_csv(infile, sep=' ', header=None)    
    except pd.errors.EmptyDataError:
        continue
    df = processing(infile)

    df.to_csv(infile.replace('.couplings', '.processed.couplings').replace('raw/', ''), index=False)

