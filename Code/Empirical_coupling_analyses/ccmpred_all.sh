#!/bin/bash

for flat_file in /stor/work/Wilke/adhock/Empirical_prot_seqs/Data/analyzed_set/*mafft.processed.flat; do
    out_file="$(basename "$flat_file" .mafft.processed.flat).ccmpred.couplings"
    echo $flat_file
    /stor/home/adhock/workspace/CCMpred/bin/ccmpred -t 50 -n 1000 $flat_file $out_file 
    wait
done
