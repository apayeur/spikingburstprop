#!/bin/bash

NUM_REAL=1
DIR="../data/noise-matching/private-noise"
mkdir -p $DIR
DATADIR="../../data/noise-matching/private-noise/"
RESULTDIR="../results/noise-matching/private-noise/"
mkdir -p $RESULTDIR

make
for W12 in 0.6
do
	for W1SOM in 0.05
	do	
		for (( i=1; i<=$NUM_REAL; i++ ))
		do
   			./sim_noise_matching_private_noise --seed $i --dir $DIR --w12 $W12 --w1som $W1SOM
		done
		
		cd ../analysis/propagation
		python plot_noisematching.py -filesuffix "_W12"$W12"_W1SOM"$W1SOM -datadir $DATADIR -resultdir "../"$RESULTDIR -numberofrealiz $NUM_REAL
		cd ../../sim
		open $RESULTDIR"Pop1_W12"$W12"_W1SOM"$W1SOM".pdf"
                open $RESULTDIR"Pop2_W12"$W12"_W1SOM"$W1SOM".pdf"
		open $RESULTDIR"SOM_W12"$W12"_W1SOM"$W1SOM".pdf"
		open $RESULTDIR"PV_W12"$W12"_W1SOM"$W1SOM".pdf"
	done
done
