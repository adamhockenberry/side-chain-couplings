#!/bin/bash

function_command='nohup python basic_simulation_include_cys.py  ../../Empirical_prot_seqs/Data/structures/1AOE_A.rosetta.pdb'

for i in $(seq 0 2); do  
    for j in $(seq 0 9); do  
        for k in $(seq -w 0 49 ); do  
            eval $function_command $i$j$k \&
        done
        wait
    
        for k in $(seq 50 99); do 
            eval $function_command $i$j$k \&
        done
        wait

    done
    wait

done
