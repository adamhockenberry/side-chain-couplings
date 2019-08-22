import pandas as pd
import glob
import shutil

db_to_assess = 'rp55'

for infile_loc in glob.glob('../Data/{}/*_{}.tsv'.format(db_to_assess, db_to_assess)):
    prot_name = infile_loc.split('/')[-1].split('_')[0] + '_' + infile_loc.split('/')[-1].split('_')[1]
    df = pd.read_csv(infile_loc, index_col=0, sep='\t', skiprows=3, skipfooter=9)
    print(prot_name, df.shape)
    if df.shape[0] >10000:
        shutil.move('../Data/fastas/{}.rosetta.fasta'.format(prot_name), '../Data/fastas/{}_only/{}.rosetta.fasta'.format(db_to_assess, prot_name))
