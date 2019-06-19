#!/bin/bash

make
current=200e-12
DIR="../data/propagation/ficurve/for-xor"
RESULT_DIR="../results/propagation/ficurve/for-xor"
FILESUFFIX=""
mkdir -p $DIR
mkdir -p $RESULT_DIR

./sim_propagation_ficurve_for_xor --somacurrent $current --dir $DIR
./sim_propagation_ficurve_for_xor --somacurrent 0. --dir $DIR
cd ../analysis/propagation/
python plot_ficurves.py  -resultdir "../"$RESULT_DIR -datadir "../"$DIR
cd ../../sim/


