#!/bin/bash

DIR="../data/xor"
mkdir -p $DIR

RESULTDIR="../results/xor"
mkdir -p $RESULTDIR

DUREX=20
ALPHA=5
NUMEX=8
NUMBERofNEURONS=200
POISSONWEIGHT=0.35

#cd ../analysis/xor
#python construct_identity_connection.py -weight $POISSONWEIGHT -numberofneuron $NUMBERofNEURONS

#cd ../../sim
make

./sim_xor --dir $DIR --durex $DUREX --alpha $ALPHA --numex $NUMEX --N $NUMBERofNEURONS
		
cd ../analysis/xor
python plot_xor.py -durex $DUREX -alpha $ALPHA -numex $NUMEX
python plot_bias.py -durex $DUREX -alpha $ALPHA -numex $NUMEX
python plot_cost.py
open ../../results/xor/Activity_Before.pdf
open ../../results/xor/Activity_DuringLearning.pdf
open ../../results/xor/Activity_After.pdf
open ../../results/xor/Wsum.pdf
cd ../../sim
