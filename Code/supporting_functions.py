import numpy as np
import pandas as pd

def process_couplings_df(df_couplings, df_contacts, primary_distance_cutoff):
    '''
    Converts the basic output of plmc into a slightly more usable dataframe
    '''
    df_couplings.columns = ['aa1_loc', 'trash1', 'aa2_loc', 'trash2', 'trash3', 'couplings']
    #First get the absolute difference in chain number between all amino acid pairs
    df_couplings['abs_diff_in_loc'] = df_couplings['aa1_loc'] - df_couplings['aa2_loc']
    df_couplings['abs_diff_in_loc'] = df_couplings['abs_diff_in_loc'].abs()
    #Only consider amino acids separated by greater than some chain distance
    df_couplings = df_couplings[df_couplings['abs_diff_in_loc'] >=primary_distance_cutoff]
    df_couplings = df_couplings.sort_values('couplings', ascending=False)

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
    df_couplings_stack = df_couplings_stack.stack().reset_index()
    df_couplings_stack.reset_index(drop=True, inplace=True)
    df_couplings_stack.columns = ['aa1_loc', 'aa2_loc', 'couplings']

    return df_couplings, df_couplings_pivot, df_couplings_stack

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
    tps = temp_df[temp_df['distance']<length_cutoff]['distance'].count()
    totals = temp_df['distance'].count()
    return  tps/totals, totals

def errors_from_df(merged_df, number_to_test, length_cutoff=8):
    merged_df.sort_values('couplings', ascending=False, inplace=True)
    temp_df = merged_df[:number_to_test]
    tps_df = temp_df[temp_df['distance']<=length_cutoff]
    tps = tps_df['distance'].count()
    fps_df = temp_df[temp_df['distance']>length_cutoff]
    fns_df = merged_df[merged_df['distance']<=length_cutoff][tps:]
    return  tps_df, fps_df, fns_df
