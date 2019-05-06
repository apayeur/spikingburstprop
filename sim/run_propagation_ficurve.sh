#!/bin/bash

make
for WSOM2 in 0.01 0.02 0.03
do
	for W12 in 0.2 0.3
	do
    		for W1SOM in 0.02 0.03 0.04
    		do
        		for current in 0 100e-12
			do
			DIR="../data/propagation/ficurve/wsom2_"$WSOM2"_w12_"$W12"_w1som_"$W1SOM
        		RESULT_DIR="../results/propagation/ficurve/wsom2_"$WSOM2
			FILESUFFIX="_w12_"$W12"_w1som_"$W1SOM
        		mkdir -p $DIR
			mkdir -p $RESULT_DIR
        		./sim_propagation_ficurve --somacurrent $current --dir $DIR --w12 $W12 --w1som $W1SOM --wsom2 $WSOM2
			done

        		cd ../analysis/propagation/
        		python plot_ficurves.py -filesuffix $FILESUFFIX -resultdir "../"$RESULT_DIR -datadir "../"$DIR
        		cd ../../sim/
    		done
	done
done


