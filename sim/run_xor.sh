#!/bin/bash

DIR="../data/xor"
mkdir -p $DIR

RESULTDIR="../results/xor"
mkdir -p $RESULTDIR

DUREX=10
ALPHA=4
NUMEX=4

make

./sim_xor --dir $DIR --durex $DUREX --alpha $ALPHA --numex $NUMEX
		
cd ../analysis/xor
python plot_xor.py -durex $DUREX -alpha $ALPHA -numex $NUMEX
open ../../results/xor/Activity.pdf
open ../../results/xor/Wsum.pdf
cd ../../sim
