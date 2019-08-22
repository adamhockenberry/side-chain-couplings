#!/bin/bash

#./workspace/hmmer-3.2.1/src/phmmer -A single_test.sto --tblout single_test.tsv -o trash.txt --incE 0.0001 -E 0.0001 single_test.fasta Protein_database_files/rp-seqs-15.fasta.gz

#./workspace/hmmer-3.2.1/easel/miniapps/esl-reformat -o single_test.afa -u afa single_test.sto

#grep -v "^#" myhits.tbl | awk '{print $1}' | esl-sfetch -f uniprot_sprot.fasta - > myhits.fa

for fasta_file in /stor/work/Wilke/adhock/Empirical_prot_seqs/Data/fastas/*.fasta; do
    prot_name=$(basename "$fasta_file" .fasta)
    echo $prot_name
    sto_name="/stor/work/Wilke/adhock/Empirical_prot_seqs/Data/rp75/${prot_name}_rp75.sto" 
    afa_name="/stor/work/Wilke/adhock/Empirical_prot_seqs/Data/rp75/${prot_name}_rp75.afa" 
    tsv_name="/stor/work/Wilke/adhock/Empirical_prot_seqs/Data/rp75/${prot_name}_rp75.tsv" 
    results_fasta_name="/stor/work/Wilke/adhock/Empirical_prot_seqs/Data/rp75/${prot_name}_rp75.fasta" 
    
    /stor/home/adhock/workspace/hmmer-3.2.1/src/phmmer -A $sto_name --tblout $tsv_name -o ./trash.txt --incE 0.0001 -E 0.0001 $fasta_file /stor/home/adhock/Protein_database_files/rp-seqs-75.fasta.gz
    wait
    
    /stor/home/adhock/workspace/hmmer-3.2.1/easel/miniapps/esl-reformat -o $afa_name -u afa $sto_name
    wait
    
    /stor/home/adhock/workspace/hmmer-3.2.1/easel/miniapps/esl-reformat -o $results_fasta_name -u fasta $sto_name
    wait
    
    rm $sto_name
    wait
done
