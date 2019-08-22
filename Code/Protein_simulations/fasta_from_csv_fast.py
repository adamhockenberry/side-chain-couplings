#!/share/apps/python-2.7.2/bin/python

import pandas as pd
from Bio import SeqIO
import glob

thresh=0.9
prot_name = '1AOE_A'
pdb_file = '/stor/work/Wilke/adhock/Empirical_prot_seqs/Data/structures/{}.rosetta.pdb'.format(prot_name)
evol_file_pattern = '/stor/work/Wilke/adhock/Protein_simulation/Results/{}/{}_thresh={}_Neff=1000_beta=1_i=*.csv'.format(prot_name, prot_name, thresh)
#nSeqs = [500, 1000]
#nSeqs = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
#nSeqs = [500, 1000, 1500, 2000, 2500, 3000]
#nSeqs = [100, 200, 300, 400]
#nMuts = [20]
nSeqs = [2000]
#nSeqs = [100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000]
#nMuts = [0.5, 1, 2, 3, 4, 5, 10]
nMuts = [10]
#nMuts = [100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
#nMuts = [100, 200]
#nMuts = [790]
#nMuts = [128, 256, 512, 1024, 2048]

    
records = list(SeqIO.parse(pdb_file, 'pdb-atom'))
assert len(records) == 1
record = records[0]
wt_seq = str(record.seq)
seq_str = str(record.seq)
nMuts = [i*len(wt_seq) for i in nMuts]

for nSeq in nSeqs:
    print(nSeq)
    seq_dicty = {}
    for i_file in glob.glob(evol_file_pattern):
        replicate_number = int(i_file.split('_i=')[-1].strip('.csv'))
        #print(i_file, replicate_number)
        if replicate_number > nSeq:
            continue
        seq_dicty[replicate_number] = {}
        seq_list = list(seq_str)
        df = pd.read_csv(i_file, index_col='Variant')
        for iteration, index in enumerate(df.index[1:nMuts[-1]+1]):
            loc = int(index[1:-1])
            assert seq_list[loc-1] == index[0],index
            seq_list[loc-1] = index[-1]
            
            if iteration+1 in nMuts:
                seq_dicty[replicate_number][iteration+1] = ''.join(seq_list)
    
    
    for nMut in nMuts:
        outfile_loc = '../Results/{}/{}_thresh={}_nMuts={}_nSeqs={}.fasta'.format(prot_name, prot_name, thresh, nMut, nSeq)
        with open(outfile_loc, 'w') as outfile:
            outfile.write('>{}\n{}\n'.format('WT', wt_seq))
            for key, val in seq_dicty.items():            
                outfile.write('>{}\n{}\n'.format(key, val[nMut]))
