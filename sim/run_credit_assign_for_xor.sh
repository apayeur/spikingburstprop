#!/bin/bash

DIR=../data/credit-assign
RESULTDIR=../results/credit-assign/for-xor
mkdir -p $DIR
mkdir -p $RESULTDIR

cd ../analysis/credit-assign
python build_current.py
cd ../../sim

make

for i in 1
	do
   		./sim_credit_assign_for_xor --seed $i --dir $DIR
	done
		
cd ../analysis/credit-assign
python plot_creditassign.py -resultdir "../"$RESULTDIR 
open ../../results/credit-assign/for-xor/CreditAssign.pdf
open ../../results/credit-assign/for-xor/Weight.pdf
cd ../../sim
