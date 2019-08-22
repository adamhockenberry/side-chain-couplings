#!/bin/bash

#for fasta_file in ../Data_new/clean_rp15_fastas/*.afa; do
#    echo "$(basename "$fasta_file" .afa).couplings" $fasta_file 
#    /stor/home/adhock/workspace/plmc/bin/plmc -n 40 -c "$(basename "$fasta_file" .afa).couplings" $fasta_file  
#    wait
#done

#for fasta_file in /stor/work/Wilke/adhock/Empirical_prot_seqs/Data/rp15/clean_fastas/*.fasta; do
#    echo "$(basename "$fasta_file" .fasta).couplings" $fasta_file 
#    /stor/home/adhock/workspace/plmc/bin/plmc -n 6 -c "$(basename "$fasta_file" .fasta).couplings" $fasta_file  
#    wait
#done

for fasta_file in /stor/work/Wilke/adhock/Protein_simulation/Results/*0.7*.fasta; do
    echo "$(basename "$fasta_file" .fasta).couplings" $fasta_file 
    /stor/home/adhock/workspace/plmc/bin/plmc -n 50 -t -1 -c "$(basename "$fasta_file" .fasta).couplings" $fasta_file  
    wait
done

#for fasta_file in ../Data/empirical_fastas/rp35/*.fasta; do
#    echo "$(basename "$fasta_file" .fasta).couplings" $fasta_file 
#    /stor/home/adhock/workspace/plmc/bin/plmc -n 40 -c "$(basename "$fasta_file" .fasta).rp35.couplings" $fasta_file  
#    wait
#done

#for fasta_file in ../Results/1B4T_A_*.fasta; do
#    echo "$(basename "$fasta_file" .fasta).couplings" $fasta_file 
#    /stor/home/adhock/workspace/plmc/bin/plmc -n 80 -c "$(basename "$fasta_file" .fasta).couplings" $fasta_file  
#    wait
#done
