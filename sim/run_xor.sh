#!/bin/bash

DIR="../data/xor"
mkdir -p $DIR

RESULTDIR="../results/xor"
mkdir -p $RESULTDIR

DUREX=20
ALPHA=5
NUMEX=2000
NUMBERofNEURONS=2000
POISSONWEIGHT=0.35

#cd ../analysis/xor
#python construct_identity_connection.py -weight $POISSONWEIGHT -numberofneuron $NUMBERofNEURONS

#cd ../../sim
make

./sim_xor --dir $DIR --durex $DUREX --alpha $ALPHA --numex $NUMEX --N $NUMBERofNEURONS
		
cd ../analysis/xor
python plot_xor_new.py -durex $DUREX -alpha $ALPHA -numex $NUMEX
open ../../results/xor/XOR_NEW.pdf
