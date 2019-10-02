#!/bin/bash

# Create directories for data and results if they don't exist already
DIR="../data/propagation"
mkdir -p $DIR
DATADIR="../../data/propagation/"
RESULTDIR="../results/propagation/"
mkdir -p $RESULTDIR

# Parameters for simulation: 
# NUM_REAL : number of realizations
# W12 : weights from population 2 to population 1 
# W1SOM : weights from SOM population to population 1
NUM_REAL=5
W12=0.37
W1SOM=0.035

# Simulate
make
for (( i=1; i<=$NUM_REAL; i++ ))
	do
   		./sim_propagation --seed $i --dir $DIR --w12 $W12 --w1som $W1SOM
	done


# Produce the figures		
cd ../analysis/propagation
python plot_propagation.py -filesuffix "_W12"$W12"_W1SOM"$W1SOM -datadir $DATADIR -resultdir "../"$RESULTDIR -numberofrealiz $NUM_REAL
python plot_raster.py -datadir $DATADIR -resultdir "../"$RESULTDIR


# Open the results
cd ../../sim
open $RESULTDIR"Pop1_W12"$W12"_W1SOM"$W1SOM".pdf"
open $RESULTDIR"Pop2_W12"$W12"_W1SOM"$W1SOM".pdf"
open $RESULTDIR"SOM_W12"$W12"_W1SOM"$W1SOM".pdf"
open $RESULTDIR"PV_W12"$W12"_W1SOM"$W1SOM".pdf"
