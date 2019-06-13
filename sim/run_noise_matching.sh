#!/bin/bash

NUM_REAL=5
DIR="../data/noise-matching"
mkdir -p $DIR
DATADIR="../../data/noise-matching/"
RESULTDIR="../results/noise-matching/"
mkdir -p $RESULTDIR

make
for W12 in 0.3
do
	for W1SOM in 0.03
	do	
		for (( i=1; i<=$NUM_REAL; i++ ))
		do
   			./sim_noise_matching --seed $i --dir $DIR --w12 $W12 --w1som $W1SOM
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
