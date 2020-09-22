# spikingburstprop
This repository contains code to simulate spiking networks with burst-mediated plasticity and multiplexing, using the Auryn simulator.

To compile this code, first compile the forked Auryn library from 
https://github.com/apayeur/auryn
and follow the instructions. 
The Auryn simulator was created by Friedeman Zenke (https://github.com/fzenke/auryn/).
The installation takes ~10 minutes on LTS Ubuntu with an Intel 5 CPU.

Python requirements:
* python 3.6.2
* numpy 1.13.1
* scipy 0.19.1
* matplotlib 2.0.2
* seaborn 0.8 
* and others

File  `./analysis/plasticity-protocol/env_py3.yml` contains the name of all packages. 
It is best to create a conda environment from the `.yml` to make sure all required Python packages are available.   

The C++ compiler is mpicxx. In the makefile `./sim/Makefile` , you will have to specify your own Auryn install path. 
This code has been tested on Mac OS Mojave 10.14.6 and Linux LTS Ubuntu.

Other, more specific, requirements are detailed below where necessary.

## Main figures
To reproduce Figs. 2-4 of our paper "*Burst-dependent synaptic plasticity can coordinate learning in hierarchical circuits*", do the following once you have cloned the present repository and with Auryn installed as described above. Runtime estimates are for a Macbook Air, 1.6 GHz Intel Core i5. 

* Figure 2B-D: From folder `./analysis/plasticity-protocol`, run the Python script `plot_all_adaptive_rule.py`.  (Runtime estimate: ~10 minutes)

* Figure 2E-G: From folder `./sim`, do: `bash run_plasticity.sh` (Runtime estimate: ~1 minute)

* Figure 3: From folder `./sim`, do: `bash run_propagation.sh` (Runtime estimate: ~1 minute)

* Figure 4: From folder `./sim`, do: `bash run_xor.sh` 
(Runtime estimate: ~3-4 hours, when `NUM_REAL=1`.) Note: You might want to change `NUM_REAL=5` to `NUM_REAL=1` in `run_xor.sh` when first testing this code. In that case, the simulation lasts 3-4 h. There is no support for parallel simulations at this point.

## Supplementary figures
Codes for the supplementatry figures are more raw than the codes for the main figures, and runtime estimates are not provided. Folder `supplementary-figures/coarse-grained-model` contains the code and data to simulate the rate model of Figure S4. Folder `supplementary-figures/supp-sim` contains the simulations for Figures S5, S10 and S11. Figures S1 and S2 can be reproduced by slightly modifying the XOR program, while Figure S3 uses `sim_xor_symfb.cpp`.

* Figure S4: From folder `supplementary-figures/coarse-grained-model/`, do: `bash run_coarse_grained.sh`. This code requires the Armadillo library from http://arma.sourceforge.net. The MNIST data are in folder `supplementary-figures/coarse-grained-model/mnist` and are provided for convenience. Code to read that data was made by Eric Yuan (see license in file `readMNIST.cpp`). 

* Figure S5: From folder `supplementary-figures/supp-sim/`, do: `bash run_dendritic_modulation.sh`.

* Figure S10: From folder `supplementary-figures/supp-sim/`, do: `bash run_normalization_of_bpcurves.sh`. Option `--w_pyr_to_pyr 0.02`, produces Figure S10b, left panel, whereas `--w_pyr_to_pyr 0.` produces Figure S10a, left panel. The right panels are produced in a similar fashion. 

* Figure S11: This one is a bit more plug-and-play. After compiling with `make`, run `./sim_direct_transmit`. Instructions to plot the results using gnuplot are provided in `sim_direct_transmit.cpp`. Several options can be tested by commenting/uncommenting relevant portions of the code.


