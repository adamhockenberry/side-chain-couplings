import requests
import pandas as pd
import glob
import time
import os.path


db = 'uniprotkb'
###First set
base_string = 'https://www.ebi.ac.uk/Tools/hmmer/download/6CB17AF2-9801-11E8-B688-A196DBC3747A.{}/score?format=fullfasta'
tsv_string = 'https://www.ebi.ac.uk/Tools/hmmer/download/6CB17AF2-9801-11E8-B688-A196DBC3747A.{}/score?format=tsv'
df = pd.read_csv('./laggard_phmmer_results.csv', index_col=0, header=None)
df.index = list(range(1, len(df.index)+1))
###Second set
#base_string = 'https://www.ebi.ac.uk/Tools/hmmer/download/926908CA-8C6F-11E8-B89C-9271E976C163.{}/score?format=fullfasta'
#df = pd.read_csv('../Data/{}_phmmer_926908CA.csv'.format(db), index_col=0, header=None)

df.columns = ['pdb', 'hit_count', 'status', 'identifier', 'description', 'e_val', 'show']

for index in df.index[:]:
#    index = str(index)
    if df.loc[index]['hit_count'] > 10000:
        continue
    prot_name = df.loc[index]['pdb'].split('.')[0]
    fname = '../../Data/{}/{}_{}.fasta.gz'.format(db, prot_name, db)
    ind_string = base_string.format(index)
    print(ind_string)
    print(index, prot_name)
    try:
        req = requests.get(ind_string, stream=True, timeout=300)
#    except requests.ConnectionError:
    except requests.exceptions.Timeout:
        print('TIMEOUT ON {}'.format(index))
        continue
    with open(fname, 'wb') as f:
        f.write(req.raw.data)
    time.sleep(30)
   
 
    ind_string = tsv_string.format(index)
    print(ind_string)
    fname = '../../Data/{}/{}_{}.tsv'.format(db, prot_name, db)
    try:
        req = requests.get(ind_string, timeout=300)
    except requests.exceptions.Timeout:
        print('TIMEOUT ON {}'.format(index))
        continue
    with open(fname, 'wb') as f:
        f.write(req.text)
    time.sleep(30)


'''
To unzip all of these files try the following commands:

gunzip *.gz

Move those files that were successful away somewhere and try the bottom code. 

It's likely that some won't work because they didnt' properly download. If this is the case, run the code below to do some clean up once or twice
'''
#
#for infile in glob.glob('../Data/{}/*.afa.gz'.format(db)):
#    prot_name = infile.split('/')[-1].split('_')[0]
#    iter_name = infile.split('/')[-1].split('_')[1]
#    print(iter_name, prot_name)
#    ind_string = base_string.format(iter_name)
#    try:
#        req = requests.get(ind_string, stream=True, timeout=300)
##    except requests.ConnectionError:
#    except requests.exceptions.Timeout:
#        print('TIMEOUT ON {}'.format(iter_name))
#        continue
#    with open('../Data/{}/{}_{}_{}.afa.gz'.format(db, prot_name, iter_name, db), 'wb') as f:
#        f.write(req.raw.data)
#
#
#
#
#
#
#
##Getting different information (not aligned fastas)
##df = pd.read_csv('../Data_new/pdb_hmmer.csv', index_col=0, header=None)
##df.columns = ['pdb', 'hit_count', 'status', 'identifier', 'description', 'e_val', 'show']
##df.head()
##
##
#base_string = 'https://www.ebi.ac.uk/Tools/hmmer/download/0384336A-3EB2-11E8-8DC8-C2F8DBC3747A.{}/score?format=json'
#db = 'pdb'
#
#for index in df.index[:]:
#    if df.loc[index]['hit_count'] < 2:
#        continue
#    ind_string = base_string.format(index)
#    print(index, df.loc[index]['pdb'])
#    try:
#        req = requests.get(ind_string, stream=True, timeout=300)
##    except requests.ConnectionError:
#    except requests.exceptions.Timeout:
#        print('TIMEOUT ON {}'.format(index))
#        continue
#    with open('../Data/{}/{}_{}_{}.json.gz'.format(db, df.loc[index]['pdb'], index, db), 'wb') as f:
#        f.write(req.raw.data)
