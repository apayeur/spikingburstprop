#!/bin/bash

DIR="../data/xor"
mkdir -p $DIR

RESULTDIR="../results/xor"
mkdir -p $RESULTDIR

make

./sim_xor --dir $DIR --durex 10 --alpha 4
		
cd ../analysis/xor
python plot_xor.py
open ../../results/xor/Activity.pdf
#open ../../results/xor/Wsum.pdf
cd ../../sim
