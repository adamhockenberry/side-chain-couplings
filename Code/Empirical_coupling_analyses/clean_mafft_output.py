from Bio import SeqIO
import glob

db = 'uniprotkb'

for original_file in glob.glob('../Data/{}/*.mafft.afa'.format(db))[:]:
    print(original_file)
    records = list(SeqIO.parse(original_file, 'fasta'))

    wt_record = records[0]

    only_seqs = []
    with open(original_file.replace('.mafft.afa', '.mafft.processed.afa'), 'w') as outfile:
        valid_indices = [i for i,j in enumerate(list(str(wt_record.seq))) if j not in ['.', '-']]
        for record in records:
            seq_list = [j for i,j in enumerate(list(str(record.seq))) if i in valid_indices]
            seq_str = ''.join(seq_list)
            seq_str = seq_str.upper()
            seq_str = seq_str.replace('.', '-').replace('X', '-')
            if seq_str not in only_seqs:
                only_seqs.append(seq_str)
                outfile.write('>{}\n{}\n'.format(record.id, seq_str))

    with open(original_file.replace('.mafft.afa', '.mafft.processed.flat'), 'w') as outfile:
        for seq in only_seqs:
            outfile.write('{}\n'.format(seq))
