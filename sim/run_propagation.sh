#!/bin/bash

# Create directories for data and results if they don't exist already
DATADIR="../data/propagation/"
mkdir -p $DATADIR
RESULTDIR="../results/propagation/"
mkdir -p $RESULTDIR

# Parameters
NUM_REAL=5  # NUM_REAL : number of realizations
W12=0.37    # W12 : weights from population 2 to population 1
W1SOM=0.035 # W1SOM : weights from SOM population to population 1

# Simulate
make
for (( i=1; i<=$NUM_REAL; i++ ))
	do
   		./sim_propagation --seed $i --dir $DATADIR --w12 $W12 --w1som $W1SOM
	done

# Produce the figures		
cd ../analysis/propagation
python plot_propagation.py -datadir "../"$DATADIR -resultdir "../"$RESULTDIR -numberofrealiz $NUM_REAL

