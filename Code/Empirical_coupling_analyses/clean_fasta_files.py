from Bio import SeqIO
import glob
from Bio.Alphabet import IUPAC

valid_aas = IUPAC.IUPACProtein.letters

db = 'uniprotkb'

for original_file in glob.glob('../Data/{}/*.fasta'.format(db))[:]:
    if '.analyze' in original_file:
        continue
    print(original_file)
    records = list(SeqIO.parse(original_file, 'fasta'))
    
    try:
        search_file = original_file.replace('{}/'.format(db), 'fastas/').replace('_{}'.format(db), '').replace('.fasta', '.rosetta.fasta')
        search_record = list(SeqIO.parse(search_file, 'fasta'))
    except IOError:
        prot_name = original_file.split('/')[-1][:6]
        other_files = glob.glob('../Data/fastas/rp15_only/*.rosetta.fasta') +\
                         glob.glob('../Data/fastas/rp35_only/*.rosetta.fasta') +\
                         glob.glob('../Data/fastas/rp55_only/*.rosetta.fasta')
        hits = [search_file for search_file in other_files if prot_name in search_file]
        assert len(hits) == 1
        search_file = hits[0]
        search_record = list(SeqIO.parse(search_file, 'fasta'))

    assert len(search_record) == 1
    records = search_record + records

    new_records = []
    all_seqs = []
    for i, record in enumerate(records):
        valid_aa_seq = True
        seq = str(record.seq).upper()
        seq_aas = list(set(seq))
        for aa in seq_aas:
            if aa not in valid_aas:
                valid_aa_seq = False
        if valid_aa_seq == False:
            continue
        if seq not in all_seqs:
            new_records.append(record)
            all_seqs.append(seq)
    

    if len(new_records) > 10000:
        new_records = new_records[:10000]
    
    with open(original_file.replace('.fasta', '.analyze.fasta'), 'w') as outfile:
        SeqIO.write(new_records, outfile, 'fasta')


#for original_file in glob.glob('../Data/rp15/*.afa')[:]:
#    print(original_file)
#    records = list(SeqIO.parse(original_file, 'fasta'))
#    wt_index = []
#    for i, record in enumerate(records):
#        if record.id.split('/')[0] in original_file:
#            wt_index.append(i)
#    assert len(wt_index) == 1
#
    #with open(original_file.replace('.afa', '.processed.afa'), 'w') as outfile:
    #    wt_record = records[wt_index[0]]
    #    valid_indices = [i for i,j in enumerate(list(str(wt_record.seq))) if j not in ['.', '-']]
    #    for record in records:
    #        seq_list = [j for i,j in enumerate(list(str(record.seq))) if i in valid_indices]
    #        outfile.write('>{}\n{}\n'.format(record.id, ''.join(seq_list)))
