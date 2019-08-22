#!/bin/bash

#/stor/home/adhock/workspace/mafft_local/bin/mafft ../Data/rp15/P37028_rp15.fasta > ./temp.out

python clean_fasta_files.py
wait

for fasta_file in /stor/work/Wilke/adhock/Empirical_prot_seqs/Data/uniprotkb/*.analyze.fasta; do
    prot_name=$(basename "$fasta_file" .analyze.fasta)
    echo $prot_name
    afa_name="/stor/work/Wilke/adhock/Empirical_prot_seqs/Data/uniprotkb/${prot_name}.mafft.afa" 
    echo $afa_name
    
    /stor/home/adhock/workspace/mafft_local/bin/mafft $fasta_file > $afa_name
    wait
done

python clean_mafft_output.py
