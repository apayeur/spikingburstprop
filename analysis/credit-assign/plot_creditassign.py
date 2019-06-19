import matplotlib.pyplot as plt
import sys
sys.path.append('../')
import raster_analysis as ra
import voltage_analysis as va
import numpy as np
import argparse
import utils_plotting as up
from build_current import example_credit_assign

parser = argparse.ArgumentParser()
parser.add_argument("-filesuffix", type=str, help="suffix indicating FF weights", default="")
args = parser.parse_args()

number_of_realizations = 1

filenames1 = []
for i in range(1, number_of_realizations+1):
    filenames1.append('../../data/credit-assign/credit_assign.0.brate1_seed'+str(i))

filenames2 = []
for i in range(1, number_of_realizations+1):
    filenames2.append('../../data/credit-assign/credit_assign.0.brate2_seed'+str(i))

filenamesPV = []
for i in range(1, number_of_realizations+1):
    filenamesPV.append('../../data/credit-assign/credit_assign.0.pvrate_seed'+str(i))

filenamesSOM = []
for i in range(1, number_of_realizations+1):
    filenamesSOM.append('../../data/credit-assign/credit_assign.0.somrate_seed'+str(i))

#filenamesVd = []
#for i in range(100):
#    filenamesVd.append('../../data/noise-matching/noise_matching{:d}.0.Vd'.format(i))

# construct currents
params = {'baseline_soma': 50.e-12,
          'baseline_dend': 75.e-12,
          'example_duration': 1.,
          'relaxation_period': 0.5,
          'number_of_examples': 3,
          'amplitudes_soma': [200.e-12, 100.e-12, -200.e-12],
          'amplitudes_dend': [200.e-12, -200.e-12, 200.e-12],  # previously +100
          'slope_duration_soma': 0. * 500.e-3,
          'slope_duration_dend': 0. * 500.e-3}
binSize = 20.e-3
times, current_soma, current_dend = example_credit_assign(params, binsize=binSize)

# output figures
inputs = {'times': times, 'dendrite': current_dend, 'soma': current_soma}
BR2 = ra.display_BRER_with_inputs(filenames2, '../../results/credit-assign/Pop2' + args.filesuffix + '.pdf', inputs,
                                  population='2', binsize=binSize, pre_stim_burn=0.)
print('size BR2 = {}'.format(len(BR2)))
print('size current dend = {}'.format(len(inputs['dendrite'])))
inputs = {'times': times, 'dendrite': BR2, 'soma': current_soma}
BR1 = ra.display_BRER_with_inputs(filenames1, '../../results/credit-assign/Pop1' + args.filesuffix + '.pdf', inputs,
                            population='1', binsize=binSize, pre_stim_burn=0.)

ra.display_rates_with_inputs(filenamesPV, '../../results/credit-assign/PV' + args.filesuffix + '.pdf', inputs, pre_stim_burn=0.)
ra.display_rates_with_inputs(filenamesSOM, '../../results/credit-assign/SOM' + args.filesuffix + '.pdf', inputs, encoding='burst', pre_stim_burn=0.)
