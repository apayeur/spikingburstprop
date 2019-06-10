#!/bin/bash

make
DIR="../data/fftransfer"
RESULT_DIR="../results/fftransfer"
mkdir -p $DIR
mkdir -p $RESULT_DIR

./sim_fftransfer --dir $DIR

cd ../analysis/fftransfer/
python plot_fftransfer.py  -resultdir "../"$RESULT_DIR -datadir "../"$DIR
open  ../../results/fftransfer/FFTransferFunction.pdf



