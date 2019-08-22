import glob
from Bio import SeqIO

re_search_records = []
for aln_file in glob.glob('../Data/analyzed_set/*.afa'):
    records = list(SeqIO.parse(aln_file, 'fasta'))
    if len(records) < 1000:
        re_search_records.append(records[0])
with open('./temp_uniprot.fasta', 'w') as outfile:
    SeqIO.write(re_search_records, outfile, 'fasta')
