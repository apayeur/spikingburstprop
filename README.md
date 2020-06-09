# spikingburstprop
This repository contains code to simulate spiking networks with burst-mediated plasticity and multiplexing, using the Auryn simulator (creator: Friedemann Zenke https://github.com/fzenke/auryn/).

To compile this code, first compile the Auryn library from 
https://github.com/apayeur/auryn
and follow the instructions.

## Main figures
To reproduce Figs. 2-4 of our paper "*Burst-dependent synaptic plasticity can coordinate learning in hierarchical circuits*", do the following once you have cloned the present repository and with Auryn installed as described above:

* Figure 2B-D: From folder `./analysis/plasticity-protocol`, run the Python script `plot_all_adaptive_rule.py`.  
This code works with the following Python packages:
	* python 3.6.2
	* numpy 1.13.1
	* scipy 0.19.1
	* matplotlib 2.0.2

	Equivalently, create a conda environment from `./analysis/plasticity-protocol/env_py3.yml`.  

* Figure 2E-G: From folder `./sim`, do: `bash run_plasticity.sh`

* Figure 3: From folder `./sim`, do: `bash run_propagation.sh`

* Figure 4: From folder `./sim`, do: `bash run_xor.sh`
(Warning!!! This can be a **long** simulation (at least on a Macbook Air, 1.6 GHz Intel Core i5). There is no support for parallel simulations at this point.)

## Supplementary figures
Codes for the supplementatry figures are more raw than the codes for the main figures. Folder `supplementary-figures/coarse-grained-model` contains the code and data to simulate the rate model of Figure S8. Folder `supplementary-figures/supp-sim` contains the simulations for Figures S1, S2 and S4. Figure S3 can be reproduced by slightly modifying the XOR program.

* Figure S1: From folder `supplementary-figures/supp-sim/`, do: `bash run_normalization_of_bpcurves.sh`. Option `--w_pyr_to_pyr 0.02`, produces Figure S1b, left panel, whereas `--w_pyr_to_pyr 0.` produces Figure S1a, left panel. The right panels are produced in a similar fashion. 

* Figure S2: This one is a bit more plug-and-play. After compiling with `make`, run `./sim_direct_transmit`. Instructions to plot the results using gnuplot are provided in `sim_direct_transmit.cpp`. Several options can be tested by commenting/uncommenting relevant portions of the code.

* Figure S4: From folder `supplementary-figures/supp-sim/`, do: `bash run_dendritic_modulation.sh`.

* Figure S8: From folder `supplementary-figures/coarse-grained-model/`, do: `bash run_coarse_grained.sh`. The MNIST data are in folder `supplementary-figures/coarse-grained-model/mnist` and are provided for convenience. Code to read that data was made by Eric Yuan (see license in file `readMNIST.cpp`). 
