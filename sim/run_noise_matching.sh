#!/bin/bash

DIR="../data/noise-matching"
mkdir -p $DIR
A=0.05

make
for W12 in 0.4
do
	for W1SOM in 0.03
	do	
		for i in 1 2 3 4 5
		do
   			./sim_noise_matching --seed $i --dir $DIR --w12 $W12 --w1som $W1SOM
		done
		
		cd ../analysis/propagation
		python plot_noisematching.py -filesuffix "_W12"$W12"_W1SOM"$W1SOM
		open ../../results/noise-matching/Pop1"_W12"$W12"_W1SOM"$W1SOM".pdf"
                open ../../results/noise-matching/Pop2"_W12"$W12"_W1SOM"$W1SOM".pdf"
		open ../../results/noise-matching/SOM"_W12"$W12"_W1SOM"$W1SOM".pdf"
		open ../../results/noise-matching/PV"_W12"$W12"_W1SOM"$W1SOM".pdf"
		cd ../../sim
	done
done
