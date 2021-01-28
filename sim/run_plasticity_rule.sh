#!/bin/bash

# Parameters
CONNECT_TYPE="AdaptiveEBCP"
TAU_PRE=100.e-3
EXAMPLE_DURATION=50
SIM_NAME="plasticity_rule"
ALPHA=10
D0=0.25 # dendritic input initial strength
W0=0.08 # initial somatic synaptic strength
NUM_REAL=1

# Create parent directory for data 
PARENT_DIR="../data/learning-rule/$CONNECT_TYPE/tau_pre$TAU_PRE"
mkdir -p $PARENT_DIR

# Create data directory for these parameter values
DATADIR="$PARENT_DIR/d0"
DATADIR+="$D0/"
DATADIR+="w0$W0"
mkdir -p $DATADIR

# Compile and simulate
make
for (( i=1; i<=$NUM_REAL; i++ ))
do
    ./sim_plasticity_rule \
    --simtime $EXAMPLE_DURATION \
    --w0 $W0 \
    --tau_pre $TAU_PRE \
    --d0 $D0 \
    --connect_type $CONNECT_TYPE \
    --dir $DATADIR \
    --sim_name $SIM_NAME \
    --alpha $ALPHA \
    --seed $i
done

# Produce the figure	
cd ../analysis/learning-rule
mkdir -p ../../results/learning-rule/Wsum/multiple-realizations
OUTFILE_NAME="../../results/learning-rule/Wsum/Wsum_vs_Time_Multiple_Realizations.pdf"
python plot_all.py -datadir "../"$DATADIR -connect_type $CONNECT_TYPE -outfile_name $OUTFILE_NAME -durex $EXAMPLE_DURATION -alpha $ALPHA -numberofrealiz $NUM_REAL
open $OUTFILE_NAME

