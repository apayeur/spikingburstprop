#!/bin/bash

# Create directories for data and results if they don't exist already
DIR="../data/xor"
mkdir -p $DIR

RESULTDIR="../results/xor"
mkdir -p $RESULTDIR

# Parameters
DUREX=20 # duration of each example
ALPHA=5 # moving average time constant
NUMEX=2000 # total number of training examples
NUMBERofNEURONS=2000 # number of neurons/population

# Compilation
make
./sim_xor --dir $DIR --durex $DUREX --alpha $ALPHA --numex $NUMEX --N $NUMBERofNEURONS

# Analysis
cd ../analysis/xor
python plot_xor_new.py -durex $DUREX -alpha $ALPHA -numex $NUMEX
open ../../results/xor/XOR_NEW.pdf
