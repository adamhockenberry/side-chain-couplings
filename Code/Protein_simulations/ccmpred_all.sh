#!/bin/bash

for flat_file in /stor/work/Wilke/adhock/Protein_simulation/Results/1AOE_A/*.flat; do
    out_file="$(basename "$flat_file" .flat).ccmpred.couplings"
    echo $flat_file
    /stor/home/adhock/workspace/CCMpred/bin/ccmpred -t 50 -n 1000 -k 10 $flat_file $out_file 
    wait
done
