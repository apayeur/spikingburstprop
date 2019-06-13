#!/bin/bash

DIR="../data/credit-assign"
mkdir -p $DIR

make

for i in 1
	do
   		./sim_credit_assign_for_xor --seed $i --dir $DIR
	done
		
cd ../analysis/credit-assign
python plot_creditassign.py
open ../../results/credit-assign/Pop1.pdf
open ../../results/credit-assign/Pop2.pdf
cd ../../sim
