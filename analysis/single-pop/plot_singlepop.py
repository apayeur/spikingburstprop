import sys
sys.path.append('../')
import raster_analysis as ra
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-filesuffix", type=str, help="suffix for file name", default="")
args = parser.parse_args()

filenames = ['../../data/single-pop/single_pop.0.brate_seed1',
             '../../data/single-pop/single_pop.0.brate_seed2',
             '../../data/single-pop/single_pop.0.brate_seed3',
             '../../data/single-pop/single_pop.0.brate_seed4',
             '../../data/single-pop/single_pop.0.brate_seed5']

# parameters
# currents in pA
max_dendritic_current = 50.
min_dendritic_current = 0.
max_somatic_current = 100.
min_somatic_current = 0.

# simulation blocks
nb_of_periods = 10
burnin = 000.e-3 # relaxation time
period = 200.e-3
segtime_maxsoma = period/2
segtime_minsoma = period/2
segtime_maxdend = 0.7*period
segtime_mindend = 0.3*period
small_overlap = 0.1*period

# construct currents
begins_dendrites = []
begins_soma = []
for i in range(10):
    begins_dendrites.append(burnin+segtime_mindend-small_overlap + i*period)
    begins_soma.append(burnin+i*period)
begins_dendrites = np.array(begins_dendrites)
begins_soma = np.array(begins_soma)

times = np.arange(1, 2., 0.001)
#current_dend = 25. + 25*np.sin(2*np.pi*times/(200e-3*np.sqrt(2.)))
current_dend = ra.alternating_square_current(times, min_dendritic_current, max_dendritic_current, begins_dendrites, segtime_maxdend)
current_soma = ra.alternating_square_current(times, min_somatic_current, max_somatic_current, begins_soma, period/2)

inputs = {'times': times, 'dendrite': current_dend, 'soma': current_soma}

# output figures
ra.display_BRER_with_inputs(filenames, '../../results/single-pop/Pop_'+args.filesuffix+'.pdf', inputs)
