#!/bin/bash

# Create directories for data and results if they don't exist already
DIR="../data/xor-test"
mkdir -p $DIR

RESULTDIR="../results/xor-test"
mkdir -p $RESULTDIR

# Parameters
DUREX=8 # duration of each example
ALPHA=2 # moving average time constant
NUMEX=2000 # total number of training examples
NUMBERofNEURONS=500 # number of neurons/population
NUM_REAL=1 # number of realizations

# Simulation
make
for (( i=1; i<=$NUM_REAL; i++ ))
    do
        ./sim_xor_test --dir $DIR --seed $i --durex $DUREX --alpha $ALPHA --numex $NUMEX --N $NUMBERofNEURONS
    done

# Analysis
#cd ../analysis/xor
#python plot_xor_new.py -durex $DUREX -alpha $ALPHA -numex $NUMEX
#open ../../results/xor/XOR_NEW.pdf
