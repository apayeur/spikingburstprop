#!/bin/bash

TAU_PRE=20.e-3
PARENT_DIR="./test-output/adaptive-ebcp"
mkdir -p $PARENT_DIR
SIM_NAME="test_adaptive_ebcp"
W0=0.06

make
function run {
		./test_adaptive_ebcp \
		--simtime 200 \
		--w0 $W0 \
        	--tau_pre $TAU_PRE \
        	--d0 $1 \
        	--dir $PARENT_DIR \
        	--sim_name $SIM_NAME
		cd test-output/adaptive-ebcp
		python plot_test_adaptive_ebcp.py -filesuffix "_d0"$1
		open Wsum_vs_Time_d0$1.pdf
		cd ../..
}


for D in 0.2
do
	run $D
done

