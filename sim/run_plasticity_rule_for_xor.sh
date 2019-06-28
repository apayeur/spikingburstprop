#!/bin/bash

CONNECT_TYPE="AdaptiveEBCP"
TAU_PRE=20.e-3
PARENT_DIR="../data/learning-rule/$CONNECT_TYPE/tau_pre$TAU_PRE"
mkdir -p ../results/learning-rule/for-xor/Wsum
mkdir -p $PARENT_DIR
EXAMPLE_DURATION=10
SIM_NAME="plasticity_rule"
W0=0.05
ALPHA=2

make
function run {
	DATADIR="$PARENT_DIR/d0"
        DATADIR+="$1/"
        DATADIR+="w0$W0"
        mkdir -p $DATADIR
	    ./sim_plasticity_rule_for_xor \
		--simtime $EXAMPLE_DURATION \
		--w0 $W0 \
        --tau_pre $TAU_PRE \
        --d0 $1 \
        --connect_type $CONNECT_TYPE \
        --dir $DATADIR \
        --sim_name $SIM_NAME \
	--alpha $ALPHA
	cd ../analysis/learning-rule
	OUTFILE_NAME="../../results/learning-rule/for-xor/Wsum/Wsum_vs_Time_"$CONNECT_TYPE"_taupre"$TAU_PRE"_w0"$W0"_alpha"$ALPHA"_durex"$EXAMPLE_DURATION".pdf"
	python plot_all.py -datadir "../"$DATADIR -connect_type $CONNECT_TYPE -outfile_name $OUTFILE_NAME -durex $EXAMPLE_DURATION -alpha $ALPHA
	open $OUTFILE_NAME
	cd ../../sim
	
}


for D in 0
do
	run $D
done

