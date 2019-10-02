import sys
sys.path.append('../')
import raster_analysis as ra
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-filesuffix", type=str, help="suffix indicating FF weights", default="")
parser.add_argument("-datadir", type=str, help="data directory", default="")
parser.add_argument("-resultdir", type=str, help="parent result directory", default="")
parser.add_argument("-numberofrealiz", type=int, help="number of realization of the noise", default=1)

args = parser.parse_args()
number_of_realizations = args.numberofrealiz

# Generate filenames to be fed to the display functions below.
filenames1 = []
for i in range(1, number_of_realizations+1):
    filenames1.append(args.datadir + 'propagation.0.brate1_seed'+str(i))

filenames2 = []
for i in range(1, number_of_realizations+1):
    filenames2.append(args.datadir + 'propagation.0.brate2_seed'+str(i))

filenamesPV = []
for i in range(1, number_of_realizations+1):
    filenamesPV.append(args.datadir + 'propagation.0.pvrate_seed'+str(i))

filenamesSOM = []
for i in range(1, number_of_realizations+1):
    filenamesSOM.append(args.datadir + 'propagation.0.somrate_seed'+str(i))


# Time and current arrays
# (amplitude and offset do not matter here because the currents are ultimately rescaled in the figure)
period = 500.e-3
period_dend = period*np.sqrt(5.)
binSize = 50.e-3
times = np.arange(0., 4., binSize)
current_soma = np.sin(2.*np.pi*times/period)
current_dend = np.sin(2.*np.pi*times/period_dend)


# Output figures
output_file_prefix = args.resultdir
_, mean, std = ra.std_BRER_over_realizations(filenames1, binsize=binSize)

inputs = {'times': times, 'dendrite': current_dend, 'soma': mean['ER']}
BR2 = ra.display_BRER_with_inputs(filenames2, output_file_prefix + 'Pop2' + args.filesuffix + '.pdf', inputs,
                            population='2', binsize=binSize)
inputs = {'times': times, 'dendrite': BR2, 'soma': current_soma}
BR1 = ra.display_BRER_with_inputs(filenames1, output_file_prefix + 'Pop1' + args.filesuffix + '.pdf', inputs,
                            population='1', binsize=binSize)

ra.display_rates_with_inputs(filenamesPV, output_file_prefix + 'PV' + args.filesuffix + '.pdf', inputs)
ra.display_rates_with_inputs(filenamesSOM, output_file_prefix + 'SOM' + args.filesuffix + '.pdf', inputs, encoding='burst')
