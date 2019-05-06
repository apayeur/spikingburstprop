import sys
sys.path.append('../analysis/')
import raster_analysis as ra
import numpy as np

# raster filenames
filenames = ['./test-output/burstpoisson.0.ras_seed1',
             './test-output/burstpoisson.0.ras_seed2',
             './test-output/burstpoisson.0.ras_seed3',
             './test-output/burstpoisson.0.ras_seed4',
             './test-output/burstpoisson.0.ras_seed5']

# parameters
# currents in pA
max_dendritic_current = 150.
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
current_dend = ra.alternating_square_current(times, min_dendritic_current, max_dendritic_current, begins_dendrites, segtime_maxdend)
current_soma = ra.alternating_square_current(times, min_somatic_current, max_somatic_current, begins_soma, period/2)

inputs = {'times': times, 'dendrite': current_dend, 'soma': current_soma}

# plot
ra.display_responses_with_inputs(filenames, './test-output/MultiplexingWithBurstPoisson.pdf', inputs)

