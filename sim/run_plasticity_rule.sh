#!/bin/bash

CONNECT_TYPE="AdaptiveEBCP"
TAU_PRE=20.e-3
PARENT_DIR="../data/learning-rule/$CONNECT_TYPE/tau_pre$TAU_PRE"
mkdir -p $PARENT_DIR
SIM_NAME="plasticity_rule"
W0=0.06

make
function run {
	DATADIR="$PARENT_DIR/d0"
        DATADIR+="$1/"
        DATADIR+="w0$W0"
        mkdir -p $DATADIR
	    ./sim_plasticity_rule \
		--simtime 20 \
		--w0 $W0 \
        --tau_pre $TAU_PRE \
        --d0 $1 \
        --connect_type $CONNECT_TYPE \
        --dir $DATADIR \
        --sim_name $SIM_NAME
        #RESULTDIR="$CONNECT_TYPE/d$1"
	cd ../analysis/learning-rule
	python plot_all.py -datadir "../"$DATADIR -connect_type $CONNECT_TYPE
	open ../../results/learning-rule/Wsum/Wsum_vs_Time_$CONNECT_TYPE.pdf
	cd ../../sim
}


for D in 0.2
do
	run $D
done

