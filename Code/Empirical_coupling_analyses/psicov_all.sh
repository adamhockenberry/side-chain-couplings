#!/bin/bash

for flat_file in /stor/work/Wilke/adhock/Empirical_prot_seqs/Data/analyzed_set/*mafft.processed.flat; do
    out_file="$(basename "$flat_file" .mafft.processed.flat).psicov.couplings"
    echo $flat_file
    /stor/home/adhock/workspace/psicov/psicov -z 50 -r 0.001 $flat_file > $out_file 
    wait
done
