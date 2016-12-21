#!/bin/bash

for N in $(seq 100 200 1500); do
	for L in $(seq 0 5 15); do
		echo $N
		echo $L
		sed -i "43s/.*/#define Nmax $N/" common.h
		sed -i "42s/.*/#define Lmax $L/" common.h
		bsub -o "DATAFILES/$N-$L.txt" 'make cleanall; make all; ./KDEMDE'
		sleep 30s #want to make sure to get the job submitted successfully first
	done
done
