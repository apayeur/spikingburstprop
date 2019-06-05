import sys
sys.path.append('../')
import raster_analysis as ra
import numpy as np

filenames1 = ['../../data/propagation/propagation.0.ras1_seed1',
              '../../data/propagation/propagation.0.ras1_seed2',
              '../../data/propagation/propagation.0.ras1_seed3',
              '../../data/propagation/propagation.0.ras1_seed4',
              '../../data/propagation/propagation.0.ras1_seed5']

filenames2 = ['../../data/propagation/propagation.0.ras2_seed1',
              '../../data/propagation/propagation.0.ras2_seed2',
              '../../data/propagation/propagation.0.ras2_seed3',
              '../../data/propagation/propagation.0.ras2_seed4',
              '../../data/propagation/propagation.0.ras2_seed5']

filenamesPV = ['../../data/propagation/propagation.0.raspv_seed1',
               '../../data/propagation/propagation.0.raspv_seed2',
               '../../data/propagation/propagation.0.raspv_seed3',
               '../../data/propagation/propagation.0.raspv_seed4',
               '../../data/propagation/propagation.0.raspv_seed5']

filenamesSOM = ['../../data/propagation/propagation.0.rassom_seed1',
               '../../data/propagation/propagation.0.rassom_seed2',
               '../../data/propagation/propagation.0.rassom_seed3',
               '../../data/propagation/propagation.0.rassom_seed4',
               '../../data/propagation/propagation.0.rassom_seed5']

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
current_dend = ra.alternating_square_current(times, min_dendritic_current, max_dendritic_current, begins_dendrites, segtime_maxdend)
current_soma = ra.alternating_square_current(times, min_somatic_current, max_somatic_current, begins_soma, period/2)

inputs = {'times': times, 'dendrite': current_dend, 'soma': current_soma}

# output figures
ra.display_BRER(filenames1, '../../results/propagation/Pop1_3examples.pdf')
ra.display_BRER(filenames2, '../../results/propagation/Pop2_3examples.pdf')
#ra.display_rates_with_inputs(filenamesPV, '../../results/propagation/PV.pdf', inputs)
#ra.display_rates_with_inputs(filenamesSOM, '../../results/propagation/SOM.pdf', inputs)
