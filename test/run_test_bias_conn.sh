#!/bin/bash

TAU_PRE=20.e-3
PARENT_DIR="./test-output/bias-conn"
mkdir -p $PARENT_DIR
SIM_NAME="test_bias_conn"
W0=0.06

make
function run {
		./test_bias_conn \
		--simtime 50 \
		--w0 $W0 \
        	--tau_pre $TAU_PRE \
        	--d0 $1 \
        	--dir $PARENT_DIR \
        	--sim_name $SIM_NAME
		cd test-output/bias-conn
		python plot_test_bias_conn.py -filesuffix "_d0"$1
		open Wsum_vs_Time_d0$1.pdf
		cd ../..
}


for D in 0.1
do
	run $D
done

