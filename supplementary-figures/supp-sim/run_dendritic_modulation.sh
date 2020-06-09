#!/bin/bash

# Create directories for data and results if they don't exist already
DIR="../supp-data/dendritic-modulation"
mkdir -p $DIR

# Compile
make

# Panel a: SOM-to-PC synaptic strength
./sim_dendritic_transfer_func_modulation --dir $DIR --wsomtopyr 0.0000 --filename "openloop"
./sim_dendritic_transfer_func_modulation --dir $DIR --wsomtopyr 0.0005 --filename "closeloop_wsomtopyr0_0005"
./sim_dendritic_transfer_func_modulation --dir $DIR --wsomtopyr 0.0010 --filename "closeloop_wsomtopyr0_001"

# Panel b: VIP-mediated disinhibition
./sim_dendritic_transfer_func_modulation --dir $DIR --vipdisinhibition 50.e-12 --filename "vipdisinhibition_50"
./sim_dendritic_transfer_func_modulation --dir $DIR --vipdisinhibition 100.e-12 --filename "vipdisinhibition_100"

# Panel c: Release probability PC-to-SOM
./sim_dendritic_transfer_func_modulation --dir $DIR --releaseprob 0.015 --filename "releaseprob0_015"
./sim_dendritic_transfer_func_modulation --dir $DIR --releaseprob 0.025 --filename "releaseprob0_025"

# Panel d: Dendritic excitability
./sim_dendritic_transfer_func_modulation --dir $DIR --burstiness 1300.e-12 --filename "burstiness1300"
./sim_dendritic_transfer_func_modulation --dir $DIR --burstiness 1400.e-12 --filename "burstiness1400"

# Plot
python ../supp-analysis/transfer.py

