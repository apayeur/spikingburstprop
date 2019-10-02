#!/bin/bash

TAU_PRE=20.e-3
PARENT_DIR="./test-output/bias-ident-conn"
mkdir -p $PARENT_DIR
SIM_NAME="test_bias_ident_conn"
W0=0.3

make
function run {
		./test_bias_ident_conn \
		--simtime 50 \
		--w0 $W0 \
        	--tau_pre $TAU_PRE \
        	--d0 $1 \
        	--dir $PARENT_DIR \
        	--sim_name $SIM_NAME
		cd test-output/bias-ident-conn
		python plot_test_bias_ident_conn.py -filesuffix "_d0"$1
		open Wsum_vs_Time_BiasIdentConn_d0$1.pdf
		cd ../..
}


for D in 0.3
do
	run $D
done

