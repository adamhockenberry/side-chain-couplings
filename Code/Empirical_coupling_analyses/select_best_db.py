import pandas as pd
import glob
import shutil

db_to_assess = 'rp55'


all_files = glob.glob('../Data/uniprotkb/*.mafft.processed.afa') +\
            glob.glob('../Data/rp75/*.mafft.processed.afa') +\
            glob.glob('../Data/rp55/*.mafft.processed.afa') +\
            glob.glob('../Data/rp35/*.mafft.processed.afa') +\
            glob.glob('../Data/rp15/*.mafft.processed.afa')

grabbed_prots = []
for ind_file in all_files:
    prot_name = ind_file.split('/')[-1].split('.')[0].split('_')
    prot_name = prot_name[0] + '_' + prot_name[1]
    print(prot_name)
    db = ind_file.split('/')[-2]
    if prot_name not in grabbed_prots:
        shutil.copyfile(ind_file, ind_file.replace('/{}/'.format(db), '/analyzed_set/'))
        flat_file = ind_file.replace('.afa', '.flat')
        shutil.copyfile(flat_file, flat_file.replace('/{}/'.format(db), '/analyzed_set/'))
        grabbed_prots.append(prot_name)            
   # prot_name = infile_loc.split('/')[-1].split('_')[0] + '_' + infile_loc.split('/')[-1].split('_')[1]
   # df = pd.read_csv(infile_loc, index_col=0, sep='\t', skiprows=3, skipfooter=9)
   # print(prot_name, df.shape)
   # if df.shape[0] >10000:
   #     shutil.copyfile('../Data/fastas/{}.rosetta.fasta'.format(prot_name), '../Data/fastas/{}_only/{}.rosetta.fasta'.format(db_to_assess, prot_name))
