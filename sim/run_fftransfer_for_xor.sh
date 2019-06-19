#!/bin/bash

make
DIR="../data/fftransfer/for-xor"
RESULT_DIR="../results/fftransfer/for-xor"
mkdir -p $DIR
mkdir -p $RESULT_DIR

./sim_fftransfer_for_xor --dir $DIR

cd ../analysis/fftransfer/
python plot_fftransfer.py  -resultdir "../"$RESULT_DIR -datadir "../"$DIR
open  ../../results/fftransfer/for-xor/FFTransferFunction.pdf



