#!/bin/bash

DIR="../data/single-pop"
mkdir -p $DIR

make
for i in 1 2 3 4 5
do
   ./sim_single_pop --seed $i --dir $DIR
done

cd ../analysis/single-pop
python plot_singlepop.py -filesuffix CONST_SOMA_CURRENT
