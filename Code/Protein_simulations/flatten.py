import glob
from Bio import SeqIO
for infile in glob.glob('../Results/1AOE_A/*.fasta'):
    print(infile)
    records = list(SeqIO.parse(infile, 'fasta'))
    with open(infile.replace('.fasta', '.flat'), 'w') as outfile:
        for record in records:
            outfile.write('{}\n'.format(str(record.seq)))

    
