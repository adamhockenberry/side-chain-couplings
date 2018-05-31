#!/usr/bin/bash


for i in $(seq 0 0); do 
    for j in $(seq 0 9); do 
        for k in $(seq -w 0 49 ); do 
            echo $i$j$k &
        done
        
        wait
        sleep 2
        
        for k in $(seq 50 99);do 
            echo $i$j$k &
        done
        
        wait
        sleep 2
    done

    wait
    sleep 2

done
