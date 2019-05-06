#!/bin/bash

DIR="../data/propagation"
mkdir -p $DIR

make
for i in 1 2 3 4 5
do
   ./sim_propagation --seed $i --dir $DIR
done

cd ../analysis/propagation
python plot_propagation.py
