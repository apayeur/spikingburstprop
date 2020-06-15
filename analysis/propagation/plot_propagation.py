import sys
sys.path.append('../')
import raster_analysis as ra
import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from utils_plotting import set_axis_color
from utils_plotting import custom_colors
from utils_plotting import remove_xticklabel
plt.style.use('../thesis_mplrc.dms')
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument("-filesuffix", type=str, help="suffix indicating FF weights", default="")
parser.add_argument("-datadir", type=str, help="data directory", default="")
parser.add_argument("-resultdir", type=str, help="parent result directory", default="")
parser.add_argument("-numberofrealiz", type=int, help="number of realization of the noise", default=1)

args = parser.parse_args()
number_of_realizations = args.numberofrealiz


# --- import data file names--- #
filenames1 = []
for i in range(1, number_of_realizations+1):
    filenames1.append(args.datadir + 'propagation.0.brate1_seed'+str(i))

filenames2 = []
for i in range(1, number_of_realizations+1):
    filenames2.append(args.datadir + 'propagation.0.brate2_seed'+str(i))

rasterfile1 = args.datadir + 'propagation.0.ras1_seed1'
rasterfile2 = args.datadir + 'propagation.0.ras2_seed1'


# --- compute mean and std ER, BR and BP --- #
binSize = 50.e-3
_, mean1, std1 = ra.std_BRER_over_realizations(filenames1, binsize=binSize)
t, mean2, std2 = ra.std_BRER_over_realizations(filenames2, binsize=binSize)
_, ER1, BR1 = ra.BRER_from_brates(args.datadir + 'propagation.0.brate1_seed1')  # for raster
_, ER2, BR2 = ra.BRER_from_brates(args.datadir + 'propagation.0.brate2_seed1')  # for raster


# --- stimulation parameters --- #
assert(abs(binSize - (t[1] - t[0])) < 1.e-6)
period = 500.e-3
period_dend = period*np.sqrt(5.)
current_soma = np.sin(2.*np.pi*t/period)
current_dend = np.sin(2.*np.pi*t/period_dend)
# Note: amplitude and offset do not matter here because the currents are ultimately rescaled in the figure)


# --- compute "input currents" for each population
pre_stim_burn = 1. # start at 1 s to avoid the transient surge of activity
min_BP1 = np.min(100*mean1['BP'][int(pre_stim_burn / binSize):]) 
max_BP1 = np.max(100*mean1['BP'][int(pre_stim_burn / binSize):])

min_BP2 = np.min(100*mean2['BP'][int(pre_stim_burn / binSize):]) 
max_BP2 = np.max(100*mean2['BP'][int(pre_stim_burn / binSize):])

min_BR2 = np.min(mean2['BR'][int(pre_stim_burn / binSize):])
max_BR2 = np.max(mean2['BR'][int(pre_stim_burn / binSize):])

min_ER1 = np.min(mean1['ER'][int(pre_stim_burn / binSize):])
max_ER1 = np.max(mean1['ER'][int(pre_stim_burn / binSize):])

min_ER2 = np.min(mean2['ER'][int(pre_stim_burn / binSize):])
max_ER2 = np.max(mean2['ER'][int(pre_stim_burn / binSize):])

min_current_dend = np.min(current_dend[int(pre_stim_burn/binSize):])
max_current_dend = np.max(current_dend[int(pre_stim_burn/binSize):])

min_current_soma = np.min(current_soma[int(pre_stim_burn/binSize):])
max_current_soma = np.max(current_soma[int(pre_stim_burn/binSize):])

# pop2 dendrite
rescaled_input_pop2_dend = min_BP2 + (current_dend - min_current_dend) * \
                                     (max_BP2 - min_BP2)/(max_current_dend - min_current_dend)

# pop2 soma
rescaled_input_pop2_soma = min_ER2 + (mean1['ER'] - min_ER1)*(max_ER2-min_ER2)/(max_ER1-min_ER1)

# pop1 dendrite
rescaled_input_pop1_dend = min_BP1 + (mean2['BR'] - min_BR2) * (max_BP1-min_BP1)/(max_BR2 - min_BR2)

# pop 1 soma
rescaled_input_pop1_soma = min_ER1 + (current_soma - min_current_soma) * \
                                     (max_ER1 - min_ER1)/(max_current_soma - min_current_soma)


# --- get event times, burst times, etc. --- #
eventtimes1, bursttimes1, _, _ = ra.get_eventtimes(rasterfile1, 0, 0.016)
eventtimes2, bursttimes2, _, _ = ra.get_eventtimes(rasterfile2, 0, 0.016)


# --- plotting --- #
output_file_prefix = args.resultdir
fig = plt.figure(figsize=(183./25.4, 183./25.4/2.))
gs = gridspec.GridSpec(4, 3, figure=fig, width_ratios=[0.6, 1. ,1.])


# empty panel on the left-hand side to accomodate the schematic
ax_empty = fig.add_subplot(gs[:, 0])
sns.despine(ax=ax_empty, bottom=True, left=True)
ax_empty.set_xticks([], [])
ax_empty.set_yticks([], [])
ax_empty.annotate('a', (-0.25*0.6, 1), xycoords='axes fraction')


# middle column: raster plots
N_neurons = 25
lineoffs = np.arange(1, N_neurons+1)
ax_raster2 = fig.add_subplot(gs[:2, 1])
ax_raster2.eventplot(eventtimes2[:N_neurons], colors=[(0, 0.45, 0.8)], lineoffsets=lineoffs, linelengths=1, lw=0.5)
ax_raster2.eventplot(bursttimes2[:N_neurons], colors=[(0.78, 0, 0)], lineoffsets=lineoffs, linelengths=1, lw=0.5)
ax_raster2.set_xlim([1, 3])
ax_raster2.set_xticks([1, 2, 3])
ax_raster2.set_yticks([1, N_neurons])
ax_raster2.spines['top'].set_visible(False)
ax_raster2.spines['right'].set_visible(False)
set_axis_color(ax_raster2, custom_colors['sky_blue'])
ax_raster2.set_ylabel('Neuron')
ax_raster2.annotate('b1', (-0.25, 1), xycoords='axes fraction')
remove_xticklabel(ax_raster2)

ax_raster2_tw = ax_raster2.twinx()
ax_raster2_tw.plot(t, 100 * BR2 / ER2, color=custom_colors['red'], lw=2)
ax_raster2_tw.plot(t, ER2, color=custom_colors['blue'], lw=2)
ax_raster2_tw.set_yticks([0, 50])
ax_raster2_tw.set_ylim([0, 60])
ax_raster2_tw.set_xlim([1, 3])
ax_raster2_tw.set_ylabel('BP [%], ER [Hz]')
remove_xticklabel(ax_raster2_tw)
set_axis_color(ax_raster2_tw, custom_colors['sky_blue'])

ax_raster1 = fig.add_subplot(gs[2:, 1])
ax_raster1.eventplot(eventtimes1[:N_neurons], colors=[(0, 0.45, 0.8)], lineoffsets=lineoffs, linelengths=1, lw=0.5)
ax_raster1.eventplot(bursttimes1[:N_neurons], colors=[(0.78, 0, 0)], lineoffsets=lineoffs, linelengths=1, lw=0.5)
ax_raster1.set_xlim([1, 3])
ax_raster1.set_xticks([1, 2, 3])
ax_raster1.set_yticks([1, N_neurons])
ax_raster1.spines['top'].set_visible(False)
ax_raster1.spines['right'].set_visible(False)
set_axis_color(ax_raster1, custom_colors['bluish_green'])
ax_raster1.set_ylabel('Neuron')
ax_raster1.set_xlabel('Time [s]')
ax_raster1.annotate('c1', (-0.25, 1), xycoords='axes fraction')

ax_raster1_tw = ax_raster1.twinx()
ax_raster1_tw.plot(t, 100 * BR1 / ER1, color=custom_colors['red'], lw=2)
ax_raster1_tw.plot(t, ER1, color=custom_colors['blue'], lw=2)
ax_raster1_tw.set_yticks([0, 50])
ax_raster1_tw.set_ylim([0, 60])
ax_raster1_tw.set_xlim([1, 3])
ax_raster1_tw.set_ylabel('BP [%], ER [Hz]')
set_axis_color(ax_raster1_tw, custom_colors['bluish_green'])


# right column: population activities
ax_pop2_dend = fig.add_subplot(gs[0, 2])
ax_pop2_dend.plot(t, 100*mean2['BP'], color=custom_colors['red'])
ax_pop2_dend.plot(t, rescaled_input_pop2_dend, '--', color=custom_colors['red'], label='$I_d$ (scaled)')
ax_pop2_dend.fill_between(t, 100*(mean2['BP'] - 2*std2['BP']), 100*(mean2['BP'] + 2*std2['BP']),
                          color=custom_colors['red'], alpha=0.5)
ax_pop2_dend.set_xlim([pre_stim_burn, 3])
ax_pop2_dend.set_ylim([0., 65])
ax_pop2_dend.set_xticks([1, 2, 3])
ax_pop2_dend.set_ylabel('BP [%]')
set_axis_color(ax_pop2_dend, custom_colors['sky_blue'])
sns.despine(ax=ax_pop2_dend)
ax_pop2_dend.legend(loc='upper right', bbox_to_anchor=(1.05, 1.3))
remove_xticklabel(ax_pop2_dend)
ax_pop2_dend.annotate('b2', (-0.25, 1.05), xycoords='axes fraction')

ax_pop2_soma = fig.add_subplot(gs[1, 2])
ax_pop2_soma.plot(t, mean2['ER'], color=custom_colors['blue'])
ax_pop2_soma.plot(t, rescaled_input_pop2_soma, '--', color=custom_colors['blue'], label='ER, pop1 (scaled)')
ax_pop2_soma.fill_between(t, mean2['ER'] - 2*std2['ER'], mean2['ER'] + 2*std2['ER'], color=custom_colors['blue'], alpha=0.5)
ax_pop2_soma.set_xlim([pre_stim_burn, 3])
ax_pop2_soma.set_ylim([0., 13.])
ax_pop2_soma.set_xticks([1, 2, 3])
ax_pop2_soma.set_ylabel('ER [Hz]')
ax_pop2_soma.legend(loc='upper right', bbox_to_anchor=(1.05, 1.3))
set_axis_color(ax_pop2_soma, custom_colors['sky_blue'])
sns.despine(ax=ax_pop2_soma)
remove_xticklabel(ax_pop2_soma)
ax_pop2_soma.annotate('b3', (-0.25, 1.05), xycoords='axes fraction')

ax_pop1_dend = fig.add_subplot(gs[2, 2])
ax_pop1_dend.plot(t, 100*mean1['BP'], color=custom_colors['red'])
ax_pop1_dend.plot(t, rescaled_input_pop1_dend, '--', color=custom_colors['orange'], label='BR, pop2 (scaled)')
ax_pop1_dend.fill_between(t, 100*(mean1['BP'] - 2*std1['BP']), 100*(mean1['BP'] + 2*std1['BP']),
                          color=custom_colors['red'], alpha=0.5)
ax_pop1_dend.set_xlim([pre_stim_burn, 3])
ax_pop1_dend.set_ylim([0., 65])
ax_pop1_dend.set_xticks([1, 2, 3])
ax_pop1_dend.set_ylabel('BP [%]')
set_axis_color(ax_pop1_dend, custom_colors['bluish_green'])
sns.despine(ax=ax_pop1_dend)
ax_pop1_dend.legend(loc='upper right', bbox_to_anchor=(1.05, 1.3))
remove_xticklabel(ax_pop1_dend)
ax_pop1_dend.annotate('c2', (-0.25, 1.05), xycoords='axes fraction')


ax_pop1_soma = fig.add_subplot(gs[3, 2])
ax_pop1_soma.plot(t, mean1['ER'], color=custom_colors['blue'])
ax_pop1_soma.plot(t, rescaled_input_pop1_soma, '--', color=custom_colors['blue'], label='$I_s$ (scaled)')
ax_pop1_soma.fill_between(t, mean1['ER'] - 2*std1['ER'], mean1['ER'] + 2*std1['ER'], color=custom_colors['blue'], alpha=0.5)
ax_pop1_soma.set_xlim([pre_stim_burn, 3])
ax_pop1_soma.set_ylim([0., 13.])
ax_pop1_soma.set_xticks([1, 2, 3])
ax_pop1_soma.set_ylabel('ER [Hz]')
ax_pop1_soma.legend(loc='upper right', bbox_to_anchor=(1.05, 1.3))
set_axis_color(ax_pop1_soma, custom_colors['bluish_green'])
sns.despine(ax=ax_pop1_soma)
ax_pop1_soma.annotate('c3', (-0.25, 1.05), xycoords='axes fraction', weight='bold')

ax_pop1_soma.set_xlabel('Time [s]')

plt.tight_layout()
plt.savefig(args.resultdir + 'fig_propagation_prelim.pdf')
plt.close()

'''
inputs = {'times': times, 'dendrite': current_dend, 'soma': mean['ER']}
BR2 = ra.display_BRER_with_inputs(filenames2, output_file_prefix + 'Pop2' + args.filesuffix + '.pdf', inputs,
                            population='2', binsize=binSize)
inputs = {'times': times, 'dendrite': BR2, 'soma': current_soma}
BR1 = ra.display_BRER_with_inputs(filenames1, output_file_prefix + 'Pop1' + args.filesuffix + '.pdf', inputs,
                            population='1', binsize=binSize)

#ra.display_rates_with_inputs(filenamesPV, output_file_prefix + 'PV' + args.filesuffix + '.pdf', inputs)
#ra.display_rates_with_inputs(filenamesSOM, output_file_prefix + 'SOM' + args.filesuffix + '.pdf', inputs, encoding='burst')
'''