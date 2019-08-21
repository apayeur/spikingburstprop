import matplotlib.pyplot as plt
import sys
sys.path.append('../')
import raster_analysis as ra
import voltage_analysis as va
import numpy as np
import argparse
import utils_plotting as up

parser = argparse.ArgumentParser()
parser.add_argument("-filesuffix", type=str, help="suffix indicating FF weights", default="")
parser.add_argument("-datadir", type=str, help="data directory", default="")
parser.add_argument("-resultdir", type=str, help="parent result directory", default="")
parser.add_argument("-numberofrealiz", type=int, help="number of realization of the noise", default=1)

args = parser.parse_args()

number_of_realizations = args.numberofrealiz

filenames1 = []
for i in range(1, number_of_realizations+1):
    filenames1.append(args.datadir + 'noise_matching.0.brate1_seed'+str(i))

filenames2 = []
for i in range(1, number_of_realizations+1):
    filenames2.append(args.datadir + 'noise_matching.0.brate2_seed'+str(i))

filenamesPV = []
for i in range(1, number_of_realizations+1):
    filenamesPV.append(args.datadir + 'noise_matching.0.pvrate_seed'+str(i))

filenamesSOM = []
for i in range(1, number_of_realizations+1):
    filenamesSOM.append(args.datadir + 'noise_matching.0.somrate_seed'+str(i))

filenamesVd = []
for i in range(100):
    filenamesVd.append(args.datadir + 'noise-matching/noise_matching{:d}.0.Vd'.format(i))

# parameters
# currents in pA
max_dendritic_current = 200.
min_dendritic_current = 100.
mean_dendritic_current = (max_dendritic_current + min_dendritic_current)/2.
A_dend = (max_dendritic_current - min_dendritic_current)/2.
max_somatic_current = 100.
min_somatic_current = 0.
mean_somatic_current = (max_somatic_current + min_somatic_current)/2.
A_soma = (max_somatic_current - min_somatic_current)/2.

# simulation blocks
nb_of_periods = 10
burnin = 000.e-3 # relaxation time
period = 500.e-3
period_dend = period*np.sqrt(5.)
segtime_maxsoma = period/2
segtime_minsoma = period/2
segtime_maxdend = 0.7*period
segtime_mindend = 0.3*period
small_overlap = 0.1*period
binSize = 50.e-3

# construct currents
begins_dendrites = []
begins_soma = []
for i in range(10):
    begins_dendrites.append(burnin+segtime_mindend-small_overlap + i*period)
    begins_soma.append(burnin+i*period)
begins_dendrites = np.array(begins_dendrites)
begins_soma = np.array(begins_soma)

times = np.arange(0, 4., binSize)
current_dend = mean_dendritic_current + A_dend*np.sin(2*np.pi*times/period_dend)
current_soma = mean_somatic_current + A_soma*np.sin(2*np.pi*times/period)

#current_dend = ra.alternating_square_current(times, min_dendritic_current, max_dendritic_current, begins_dendrites, segtime_maxdend)
#current_soma = ra.alternating_square_current(times, min_somatic_current, max_somatic_current, begins_soma, period/2)


'''
# display voltage versus time
t, Vd = va.population_voltage(filenamesVd)
t2, ER, BR = ra.BRER_from_brates(filenames2[0])
plt.figure(figsize=(6, 2.5))
plt.subplot(121)
plt.plot(t[t>1.], up.adjust_range(Vd[t>1.], 100*BR[t2>1]/ER[t2>1]), label=r'$\bar{V}_d$ (rescaled)')
plt.plot(t2[t2>1], 100*BR[t2>1]/ER[t2>1], label='BP [\%]')
plt.plot(t2[t2>1], up.adjust_range(BR[t2>1], 100*BR[t2>1]/ER[t2>1]), label='BR (rescaled)')
plt.xlabel('Time [s]')
plt.legend()
#plt.ylabel(r'$\bar{V}_d$ (rescaled)')

plt.subplot(122)
plt.plot(1000*Vd[::int(20./0.1)], 100*BR/ER, 'o', markerfacecolor='white', markeredgewidth=2, alpha=0.5)
plt.xlabel(r'$\bar{V}_d$ [mV]')
plt.ylabel('BP [\%]')
plt.tight_layout()
plt.show()
'''
# output figures
output_file_prefix = args.resultdir

inputs = {'times': times, 'dendrite': current_dend, 'soma': current_soma}
BR2 = ra.display_BRER_with_inputs(filenames2, output_file_prefix + 'Pop2' + args.filesuffix + '.pdf', inputs,
                            population='2', binsize=binSize)
inputs = {'times': times, 'dendrite': BR2, 'soma': current_soma}
BR1 = ra.display_BRER_with_inputs(filenames1, output_file_prefix + 'Pop1' + args.filesuffix + '.pdf', inputs,
                            population='1', binsize=binSize)

ra.display_rates_with_inputs(filenamesPV, output_file_prefix + 'PV' + args.filesuffix + '.pdf', inputs)
ra.display_rates_with_inputs(filenamesSOM, output_file_prefix + 'SOM' + args.filesuffix + '.pdf', inputs, encoding='burst')
