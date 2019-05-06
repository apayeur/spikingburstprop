import sys
sys.path.append('../')
import raster_analysis as ra
import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-filesuffix", type=str, help="suffix indicating FF weights", default="")
args = parser.parse_args()

filenames1 = ['../../data/noise-matching/noise_matching.0.brate1_seed1']
              #'../../data/noise-matching/noise_matching.0.brate1_seed2',
              #'../../data/noise-matching/noise_matching.0.brate1_seed3',
              #'../../data/noise-matching/noise_matching.0.brate1_seed4',
              #'../../data/noise-matching/noise_matching.0.brate1_seed5']

filenames2 = ['../../data/noise-matching/noise_matching.0.brate2_seed1']
              #'../../data/noise-matching/noise_matching.0.brate2_seed2',
              #'../../data/noise-matching/noise_matching.0.brate2_seed3',
              #'../../data/noise-matching/noise_matching.0.brate2_seed4',
              #'../../data/noise-matching/noise_matching.0.brate2_seed5']

filenamesPV = ['../../data/noise-matching/noise_matching.0.pvrate_seed1']
               #'../../data/noise-matching/noise_matching.0.pvrate_seed2',
               #'../../data/noise-matching/noise_matching.0.pvrate_seed3',
               #'../../data/noise-matching/noise_matching.0.pvrate_seed4',
               #'../../data/noise-matching/noise_matching.0.pvrate_seed5']

filenamesSOM = ['../../data/noise-matching/noise_matching.0.somrate_seed1']
               #'../../data/noise-matching/noise_matching.0.somrate_seed2',
               #'../../data/noise-matching/noise_matching.0.somrate_seed3',
               #'../../data/noise-matching/noise_matching.0.somrate_seed4',
               #'../../data/noise-matching/noise_matching.0.somrate_seed5']

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
ra.display_BRER_with_inputs(filenames1, '../../results/noise-matching/Pop1' + args.filesuffix + '.pdf', inputs)
ra.display_BRER_with_inputs(filenames2, '../../results/noise-matching/Pop2' + args.filesuffix + '.pdf', inputs)
ra.display_rates_with_inputs(filenamesPV, '../../results/noise-matching/PV' + args.filesuffix + '.pdf', inputs)
ra.display_rates_with_inputs(filenamesSOM, '../../results/noise-matching/SOM' + args.filesuffix + '.pdf', inputs)
