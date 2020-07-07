#!/bin/bash

# Create directories for data and results if they don't exist already
DIR="../data/propagation"
mkdir -p $DIR
DATADIR="../../data/propagation/"
RESULTDIR="../results/propagation/inverted/"
mkdir -p $RESULTDIR

# Parameters
NUM_REAL=5 # NUM_REAL : number of realizations
W12=0.37 # W12 : weights from population 2 to population 1
W1SOM=0.035 # W1SOM : weights from SOM population to population 1

# Simulate
make
for (( i=1; i<=$NUM_REAL; i++ ))
	do
   		./sim_propagation_inverted --seed $i --dir $DIR --w12 $W12 --w1som $W1SOM
	done

# Produce the figures		
cd ../analysis/propagation
python plot_propagation_inverted.py -filesuffix "_W12"$W12"_W1SOM"$W1SOM -datadir $DATADIR -resultdir "../"$RESULTDIR -numberofrealiz $NUM_REAL

