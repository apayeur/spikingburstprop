#!/bin/bash

CONNECT_TYPE="EBCP"
TAU_PRE=20.e-3
PARENT_DIR="../data/learning-rule/$CONNECT_TYPE/tau_pre$TAU_PRE"
mkdir -p $PARENT_DIR
SIM_NAME="plasticity_rule"
W0=0.03

make
function run {

	    ./sim_plasticity_rule \
		--simtime 200 \
		--w0 $W0 \
        --tau_pre $TAU_PRE \
        --d0 $1 \
        --connect_type $CONNECT_TYPE \
        --dir "." \
        --sim_name $SIM_NAME
        #RESULTDIR="$CONNECT_TYPE/d$1"
        DATADIR="$PARENT_DIR/d0"
        DATADIR+="$1/"
	DATADIR+="w0$W0"
        mkdir -p $DATADIR
        cp $SIM_NAME* $DATADIR
        rm $SIM_NAME*
}


for D in 0.2
do
	run $D
done

