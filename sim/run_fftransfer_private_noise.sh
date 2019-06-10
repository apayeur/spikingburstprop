#!/bin/bash

make
DIR="../data/fftransfer/private-noise"
RESULT_DIR="../results/fftransfer/private-noise"
mkdir -p $DIR
mkdir -p $RESULT_DIR

./sim_fftransfer_private_noise --dir $DIR

cd ../analysis/fftransfer/
python plot_fftransfer.py  -resultdir "../"$RESULT_DIR -datadir "../"$DIR
open  ../../results/fftransfer/private-noise/FFTransferFunction.pdf



