#!/bin/bash

for c in $(seq 0.1 0.02 0.31); do
        while read k; do
                bsub -q week -o "DATAFILES/$c-$k.txt" ./KDEMDE $k $c
        done <k_vals.txt
done
