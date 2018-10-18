import numpy as np
import pandas as pd

def process_couplings_df(df_couplings, df_contacts, primary_distance_cutoff):
    '''
    Converts the basic output of plmc into a slightly more usable dataframe
    '''
    df_couplings_pivot = df_couplings.pivot(index='aa1_loc', columns='aa2_loc', values='couplings')    
    
    for col in df_contacts.columns:
        if col not in df_couplings_pivot.columns:
            df_couplings_pivot[col] = np.nan
    for index in df_contacts.index:
        if index not in df_couplings_pivot.index:
             df_couplings_pivot = df_couplings_pivot.append(pd.Series(name=index))
    
    df_couplings_pivot.sort_index(inplace=True)
    df_couplings_pivot = df_couplings_pivot.reindex_axis(sorted(df_couplings_pivot.columns), axis=1)
        
    df_couplings_stack = df_couplings_pivot.where(np.triu(np.ones(df_couplings_pivot.shape)).astype(np.bool))
    df_couplings_stack = df_couplings_stack.stack(dropna=False).reset_index()
    df_couplings_stack = df_couplings_stack[df_couplings_stack['aa1_loc'] < df_couplings_stack['aa2_loc']]
    df_couplings_stack.reset_index(drop=True, inplace=True)
    df_couplings_stack.columns = ['aa1_loc', 'aa2_loc', 'couplings']
    
    df_couplings_stack['abs_diff_in_loc'] = df_couplings_stack['aa1_loc'] - df_couplings_stack['aa2_loc']
    df_couplings_stack['abs_diff_in_loc'] = df_couplings_stack['abs_diff_in_loc'].abs()
    ##Only consider amino acids separated by greater than some chain distance
    df_couplings_stack = df_couplings_stack[df_couplings_stack['abs_diff_in_loc'] >=primary_distance_cutoff]
    df_couplings_stack.reset_index(drop=True, inplace=True)

    return df_couplings_stack, df_couplings_pivot

def process_contacts_df(df_contacts, primary_distance_cutoff):
    df_contacts.columns = df_contacts.columns.astype(int)
    df_contacts_stack = df_contacts.where(np.triu(np.ones(df_contacts.shape)).astype(np.bool))
    df_contacts_stack = df_contacts_stack.stack().reset_index()

    df_contacts_stack['abs_diff'] = df_contacts_stack['level_0'] - df_contacts_stack['level_1']
    df_contacts_stack['abs_diff'] = df_contacts_stack['abs_diff'].abs()
    df_contacts_stack = df_contacts_stack[df_contacts_stack['abs_diff'] >= primary_distance_cutoff]
    df_contacts_stack.reset_index(drop=True, inplace=True)
    df_contacts_stack.columns = ['aa1_loc', 'aa2_loc', 'distance', 'primary_chain_distance']
    return df_contacts, df_contacts_stack

def process_angles_df(df_angles, primary_distance_cutoff):
    primary_distance_cutoff = 12

    df_angles.fillna(value=-100, inplace=True)
    df_angles.columns = df_angles.columns.astype(int)
    df_angles_stack = df_angles.where(np.triu(np.ones(df_angles.shape)).astype(np.bool))
    df_angles_stack = df_angles_stack.stack().reset_index()
    
    
    df_angles_stack['abs_diff'] = df_angles_stack['level_0'] - df_angles_stack['level_1']
    df_angles_stack['abs_diff'] = df_angles_stack['abs_diff'].abs()
    df_angles_stack = df_angles_stack[df_angles_stack['abs_diff'] >= primary_distance_cutoff]
    df_angles_stack.reset_index(drop=True, inplace=True)
    df_angles_stack.columns = ['aa1_loc', 'aa2_loc', 'angles1', 'primary_chain_distance']
    
    df_angles2_stack = df_angles.where(np.tril(np.ones(df_angles.shape)).astype(np.bool))
    df_angles2_stack = df_angles2_stack.stack().reset_index()
    
    
    df_angles2_stack['abs_diff'] = df_angles2_stack['level_0'] - df_angles2_stack['level_1']
    df_angles2_stack['abs_diff'] = df_angles2_stack['abs_diff'].abs()
    df_angles2_stack = df_angles2_stack[df_angles2_stack['abs_diff'] >= primary_distance_cutoff]
    df_angles2_stack.reset_index(drop=True, inplace=True)
    df_angles2_stack.columns = ['aa1_loc', 'aa2_loc', 'angles2', 'primary_chain_distance']
    df_angles2_stack = df_angles2_stack.sort_values(['aa2_loc', 'aa1_loc'])
    df_angles2_stack.rename(dict(zip(list(df_angles2_stack.index), list(df_angles_stack.index))), inplace=True)
    
    assert all(df_angles_stack['aa1_loc'] == df_angles2_stack['aa2_loc'])
    assert all(df_angles_stack['aa2_loc'] == df_angles2_stack['aa1_loc'])
    df_angles_stack = pd.concat([df_angles_stack, df_angles2_stack[['angles2']]],\
                                axis=1, join_axes=[df_angles_stack.index])
    df_angles.replace(-100., np.nan, inplace=True)
    df_angles_stack['angles1'] = df_angles_stack['angles1'].replace(-100, np.nan)
    df_angles_stack['angles2'] = df_angles_stack['angles2'].replace(-100, np.nan)
    return df_angles_stack

def merge_contacts_couplings(df_contacts_stack, df_couplings_stack, seq):
    assert all(df_couplings_stack['aa1_loc'] == df_contacts_stack['aa1_loc'])
    assert all(df_couplings_stack['aa2_loc'] == df_contacts_stack['aa2_loc'])
    merged_df = pd.concat([df_contacts_stack, df_couplings_stack[['couplings']]],\
                          axis=1, join_axes=[df_contacts_stack.index])
    merged_df.sort_values('couplings', ascending=False, inplace=True)
    
    seq_dict = dict([(i+1,j) for i,j in list(enumerate(list(seq)))])
    merged_df['aa1_aa'] = merged_df['aa1_loc'].map(seq_dict)
    merged_df['aa2_aa'] = merged_df['aa2_loc'].map(seq_dict)
    
    return merged_df
    
def ppv_from_df(merged_df, number_to_test, length_cutoff=8):
    temp_df = merged_df[:number_to_test]
    tps = temp_df[temp_df['distance']<=length_cutoff]['distance'].count()
    totals = temp_df['distance'].count()
    return  tps/totals, totals

def recall_from_df(merged_df, number_to_test, length_cutoff=8):
    temp_df = merged_df[:number_to_test]
    tps = temp_df[temp_df['distance']<=length_cutoff]['distance'].count()
    temp_df = merged_df[number_to_test:]
    fns = temp_df[temp_df['distance']<=length_cutoff]['distance'].count()
    totals = tps+fns
    return  tps/totals, totals 
