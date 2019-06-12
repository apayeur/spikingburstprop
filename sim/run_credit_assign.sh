#!/bin/bash

DIR="../data/credit-assign"
mkdir -p $DIR
A=0.05

make
for W12 in 0.4
do
	for W1SOM in 0.03
	do	
		for i in 1
		do
   			./sim_credit_assign --seed $i --dir $DIR --w12 $W12 --w1som $W1SOM
		done
		
		cd ../analysis/credit-assign
		python plot_creditassign.py -filesuffix "_W12"$W12"_W1SOM"$W1SOM
		open ../../results/credit-assign/Pop1"_W12"$W12"_W1SOM"$W1SOM".pdf"
                open ../../results/credit-assign/Pop2"_W12"$W12"_W1SOM"$W1SOM".pdf"
		cd ../../sim
	done
done
