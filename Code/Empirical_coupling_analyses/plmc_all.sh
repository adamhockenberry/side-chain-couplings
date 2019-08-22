#!/bin/bash

for afa_file in /stor/work/Wilke/adhock/Empirical_prot_seqs/Data/analyzed_set/*mafft.processed.afa; do
    out_file="$(basename "$afa_file" .mafft.processed.afa).plmc.couplings"
    echo $afa_file
    /stor/home/adhock/workspace/plmc/bin/plmc -n 50 -c $out_file $afa_file  
    wait
done
