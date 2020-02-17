# spikingburstprop
This repository contains code to simulate spiking networks with burst-mediated plasticity and multiplexing, using the Auryn simulator (creator: Friedemann Zenke https://github.com/fzenke/auryn/).

To compile this code, first compile the auryn library from 
https://github.com/apayeur/auryn
and follow the instructions.

To reproduce Figs. 2-3-4 of our paper "*Burst-dependent synaptic plasticity can coordinate learning in hierarchical circuits*", do the following once you have cloned the present repositiory and with Auryn installed as described above:

* Figure 2B-D:
From folder `./analysis/plasticity-protocol`, run the Python script `plot_all_adaptive_rule.py`. This code produces a `.pdf` in folder `../../results/learning-rule`, so make sure that the latter exists before running the code.  
This code works with the following Python packages:
** python 3.6.2
** numpy 1.13.1
** scipy 0.19.1
** matplotlib 2.0.2
Equivalently, create a conda environment from `./analysis/plasticity-protocol/env_py3.yml`.  

* Figure 2E-G:
From folder `./sim`, do: `bash run_plasticity.sh`

* Figure 3:
From folder `./sim`, do: `bash run_propagation.sh`

* Figure 4:
From forder `./sim`, do: `bash run_xor.sh`
(Warning!!! This is a **long** simulation (48h + on a Macbook Air, 1.6 GHz Intel Core i5). There is no support for parallel simulations at this point.)



