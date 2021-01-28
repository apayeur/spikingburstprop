import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import sys
sys.path.append('../')
plt.style.use('../thesis_mplrc.dms')
from utils_plotting import remove_xticklabel
from utils_plotting import custom_colors
from utils_plotting import set_axis_color
import raster_analysis as ra
#from raster_analysis import add_step
import seaborn as sns

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-durex", type=float, help="duration of example")
parser.add_argument("-alpha", type=float, help="moving average time scale")
parser.add_argument("-numex", type=int, help="number of training examples")
parser.add_argument("-numberofrealiz", type=int, help="number of realization of the noise", default=1)

args = parser.parse_args()
number_of_realizations = args.numberofrealiz


# import data
data_file_prefix = '../../data/xor-test/sym-ff/multiple-reals/xor-test.0.brate_'
datadir = '../../data/xor-test/sym-ff/multiple-reals/'
filenames_output = []
for i in range(1, number_of_realizations+1):
    filenames_output.append(data_file_prefix + 'output_seed'+str(i))
output_rates = np.loadtxt(data_file_prefix + 'output_seed1')
hidden1_rates = np.loadtxt(data_file_prefix + 'hidden1_seed1')
hidden2_rates = np.loadtxt(data_file_prefix + 'hidden2_seed1')
input1_rates = np.loadtxt(data_file_prefix + 'input1_seed1')
input2_rates = np.loadtxt(data_file_prefix + 'input2_seed1')
#cost_epoch = np.loadtxt(datadir + 'cost_epoch.dat')


# compute output rate
binSize = hidden1_rates[1,0] - hidden1_rates[0,0]
time_output, mean_output, std_output = ra.std_BRER_over_realizations(filenames_output, binsize=binSize)


# define selected pre, during and post learning intervals
after_learning = [hidden1_rates.shape[0]-int(4*args.durex/binSize), hidden1_rates.shape[0]+1]
before_learning = [int(3*args.alpha/binSize)+1, int(3*args.alpha/binSize)+int(4*args.durex/binSize)+1]
during_learning_interval = [3*args.alpha + args.numex//5*args.durex,
                            3*args.alpha + args.numex//5*args.durex + 4*args.durex]
during_learning = [int(during_learning_interval[0]/binSize) + 1,
                   int(during_learning_interval[1]/binSize) + 1]

er_trace = np.loadtxt(datadir + 'xor-test0.0.trevent_seed1')
br_trace = np.loadtxt(datadir + 'xor-test0.0.trburst_seed1')
for i in range(1,50):
    er_tmp = np.loadtxt(datadir + 'xor-test'+str(i)+'.0.trevent_seed1')
    br_tmp = np.loadtxt(datadir + 'xor-test'+str(i)+'.0.trburst_seed1')
    er_trace[:,1] += er_tmp[:,1]
    br_trace[:,1] += br_tmp[:,1]
er_trace[:,1] /= 50.
br_trace[:,1] /= 50.
er_trace = er_trace
br_trace = br_trace

e_min = 2.
e_max = 10.
e_th = e_min + 0.5*(e_max - e_min)



fig = plt.figure(figsize=(183/25.4, 4.5/7.5*183/25.4))
# first subplot is empty

# inputs
ax1 = fig.add_subplot(337)
ax1.plot(input1_rates[before_learning[0]:before_learning[1], 0],
         input1_rates[before_learning[0]:before_learning[1], 2],
         color=custom_colors['blue'], label='input 1')
ax1.plot(input2_rates[before_learning[0]:before_learning[1], 0],
         input2_rates[before_learning[0]:before_learning[1], 2], '-', lw=1, color=custom_colors['blue'], label='input 2')
xtik = list(3*args.alpha + args.durex*np.arange(5))
ax1.text(xtik[0] + args.durex/2., -0.5, '(0, 0)', horizontalalignment='center', verticalalignment='top')
ax1.text(xtik[1] + args.durex/2., -0.5, '(1, 0)', horizontalalignment='center', verticalalignment='top')
ax1.text(xtik[2] + args.durex/2., -0.5, '(0, 1)', horizontalalignment='center', verticalalignment='top')
ax1.text(xtik[3] + args.durex/2., -0.5, '(1, 1)', horizontalalignment='center', verticalalignment='top')
ax1.set_ylabel('ER, input [Hz]')
ax1.set_ylim([0, 11])
#ax1.set_title('Inputs for all examples', fontdict=None, loc='center', fontsize=9)
ax1.set_xticks(xtik)
ax1.set_yticks([0, 10])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.legend(loc='best')
remove_xticklabel(ax1)


# output before/after
ax2 = fig.add_subplot(332)
ax2.plot(time_output[before_learning[0]:before_learning[1]],
         mean_output['ER'][before_learning[0]:before_learning[1]],
         linestyle='dashed', color=custom_colors['blue'], label='before')
ax2.fill_between(time_output[before_learning[0]:before_learning[1]],
                 mean_output['ER'][before_learning[0]:before_learning[1]] - 2*std_output['ER'][before_learning[0]:before_learning[1]],
                 mean_output['ER'][before_learning[0]:before_learning[1]] + 2*std_output['ER'][before_learning[0]:before_learning[1]],
                 color=custom_colors['blue'], alpha=0.5, lw=0)
ax2.plot(time_output[before_learning[0]:before_learning[1]],
         mean_output['ER'][after_learning[0]:after_learning[1]],
         color=custom_colors['blue'], label='after')
ax2.fill_between(time_output[before_learning[0]:before_learning[1]],
                 mean_output['ER'][after_learning[0]:after_learning[1]] - 2*std_output['ER'][after_learning[0]:after_learning[1]],
                 mean_output['ER'][after_learning[0]:after_learning[1]] + 2*std_output['ER'][after_learning[0]:after_learning[1]],
                 color=custom_colors['blue'], alpha=0.5, lw=0)
ax2.plot(time_output[before_learning[0]:before_learning[1]],
         e_th*np.ones(np.shape(mean_output['ER'][before_learning[0]:before_learning[1]])), 'k:', lw=0.5)
ax2.text(xtik[0] + args.durex/2., -0.5, '(0, 0)', horizontalalignment='center', verticalalignment='top')
ax2.text(xtik[1] + args.durex/2., -0.5, '(1, 0)', horizontalalignment='center', verticalalignment='top')
ax2.text(xtik[2] + args.durex/2., -0.5, '(0, 1)', horizontalalignment='center', verticalalignment='top')
ax2.text(xtik[3] + args.durex/2., -0.5, '(1, 1)', horizontalalignment='center', verticalalignment='top')
ax2.set_ylabel('ER, output [Hz]')
ax2.set_ylim([0, 11])
#ax2.set_title('Output before/after learning', fontdict=None, loc='center', fontsize=9)
ax2.set_xticks(xtik)
ax2.set_yticks([0, 10])
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.legend(loc='best')
remove_xticklabel(ax2)
#set_axis_color(ax2, custom_colors['sky_blue'])


# error generation at output layer
ax3 = fig.add_subplot(333)
r_out = output_rates[during_learning[0]:during_learning[1], :]
teacher = np.zeros(np.shape(r_out[:,2]))
er_est = er_trace[during_learning[0]:during_learning[1], 1]/args.alpha
xtik = list(during_learning_interval[0] + args.durex*np.arange(5))

ra.add_step(teacher, int(0.9*args.durex/binSize), int(0.1*args.durex/binSize)-1, -1./(e_max - er_est[int(0.9*args.durex/binSize)]))
ra.add_step(teacher, int(1.9*args.durex/binSize), int(0.1*args.durex/binSize)-1, 1./(er_est[int(1.9*args.durex/binSize)] - e_min))
ra.add_step(teacher, int(2.9*args.durex/binSize), int(0.1*args.durex/binSize)-1, 1./(er_est[int(2.9*args.durex/binSize)] - e_min))
ra.add_step(teacher, int(3.9*args.durex/binSize), int(0.1*args.durex/binSize)-1, -1./(e_max - er_est[int(3.9*args.durex/binSize)]))

ax3.plot(r_out[:, 0], 100*r_out[:, 1]/r_out[:, 2], color=custom_colors['red'])
ax3.text(xtik[0] + args.durex/2., -4.5, '(0, 0)', horizontalalignment='center', verticalalignment='top')
ax3.text(xtik[1] + args.durex/2., -4.5, '(1, 0)', horizontalalignment='center', verticalalignment='top')
ax3.text(xtik[2] + args.durex/2., -4.5, '(0, 1)', horizontalalignment='center', verticalalignment='top')
ax3.text(xtik[3] + args.durex/2., -4.5, '(1, 1)', horizontalalignment='center', verticalalignment='top')
#ax3.set_title('Error generation at output layer', fontdict=None, loc='center', fontsize=9)
remove_xticklabel(ax3)
#set_axis_color(ax3, custom_colors['sky_blue'])
ax3.set_yticks([0, 45, 90])
ax3.set_ylabel('BP, output [%]', color=custom_colors['red'])
ax3.spines['top'].set_visible(False)

ax3_tw = ax3.twinx()
ax3_tw.plot(r_out[:, 0], teacher, '--', lw=1, color=custom_colors['pink'])
ax3_tw.set_xticks(xtik)
#set_axis_color(ax3_tw, custom_colors['sky_blue'])
#ax3_tw.set_ylim(ax3_tw.get_yticks()[0], ax3_tw.get_yticks()[-1])
ax3_tw.set_yticks([-0.4, 0, 0.4])
ax3_tw.set_ylabel('Teaching current [nA]', color=custom_colors['pink'])
ax3_tw.spines['top'].set_visible(False)


# error propagation in hidden layers
ax4 = fig.add_subplot(335)
r_hid1 = hidden1_rates[during_learning[0]:during_learning[1], :]
r_out = output_rates[during_learning[0]:during_learning[1], :]
BP = 100*r_hid1[:, 1]/r_hid1[:, 2]
min_BP = np.min(BP)
max_BP = np.max(BP)
min_current_dend = np.min(r_out[:,1])
max_current_dend = np.max(r_out[:,1])
rescaled_BR = min_BP + (r_out[:,1] - min_current_dend)*(max_BP-min_BP)/(max_current_dend-min_current_dend)
ax4.plot(r_hid1[:, 0], BP, color=custom_colors['red'])
ax4.plot(r_out[:, 0], rescaled_BR, ':', color=custom_colors['orange'], label='output BR (rescaled)')
ax4.set_xticks(xtik)
ax4.set_yticks([0, 50])
ax4.set_ylim([0, 60])
ax4.set_ylabel('BP, hidden 1 [%]')
ax4.text(xtik[0] + args.durex/2., -2.5, '(0, 0)', horizontalalignment='center', verticalalignment='top')
ax4.text(xtik[1] + args.durex/2., -2.5, '(1, 0)', horizontalalignment='center', verticalalignment='top')
ax4.text(xtik[2] + args.durex/2., -2.5, '(0, 1)', horizontalalignment='center', verticalalignment='top')
ax4.text(xtik[3] + args.durex/2., -2.5, '(1, 1)', horizontalalignment='center', verticalalignment='top')
#ax4.set_title('Error propagation to hidden 1', fontdict=None, loc='center', fontsize=9)
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.legend(loc='best')
remove_xticklabel(ax4)
#set_axis_color(ax4, custom_colors['bluish_green'])


ax5 = fig.add_subplot(336)
r_hid2 = hidden2_rates[during_learning[0]:during_learning[1], :]
BP = 100*r_hid2[:, 1]/r_hid2[:, 2]
min_BP = np.min(BP)
max_BP = np.max(BP)
BR_out = -output_rates[during_learning[0]:during_learning[1], 1]
min_current_dend = np.min(BR_out)
max_current_dend = np.max(BR_out)
rescaled_BR = min_BP + (BR_out - min_current_dend)*(max_BP-min_BP)/(max_current_dend-min_current_dend)
ax5.plot(r_hid2[:, 0], BP, color=custom_colors['red'])
ax5.plot(r_out[:, 0], rescaled_BR, ':', color=custom_colors['orange'], label='output BR (inverted & rescaled)')

ax5.set_xticks(list(during_learning[0]*binSize + args.durex*np.arange(5)-1))
ax5.set_yticks([0, 50])
ax5.set_ylim([0, 60])
ax5.set_ylabel('BP, hidden 2 [%]')
ax5.text(xtik[0] + args.durex/2., -2.5, '(0, 0)', horizontalalignment='center', verticalalignment='top')
ax5.text(xtik[1] + args.durex/2., -2.5, '(1, 0)', horizontalalignment='center', verticalalignment='top')
ax5.text(xtik[2] + args.durex/2., -2.5, '(0, 1)', horizontalalignment='center', verticalalignment='top')
ax5.text(xtik[3] + args.durex/2., -2.5, '(1, 1)', horizontalalignment='center', verticalalignment='top')
#ax5.set_title('Error propagation to hidden 2', fontdict=None, loc='center', fontsize=9)
ax5.spines['top'].set_visible(False)
ax5.spines['right'].set_visible(False)
ax5.legend(loc='lower left')
remove_xticklabel(ax5)
#set_axis_color(ax5, custom_colors['bluish_green'])

# hidden 1 before/after
ax6 = fig.add_subplot(338)
ax6.plot(hidden1_rates[before_learning[0]:before_learning[1], 0], hidden1_rates[before_learning[0]:before_learning[1], 2],
         linestyle='dashed', color=custom_colors['blue'], label='before')
ax6.plot(hidden1_rates[before_learning[0]+1:before_learning[1], 0], hidden1_rates[after_learning[0]+1:after_learning[1], 2],
         color=custom_colors['blue'], label='after')
xtik = list(3*args.alpha + args.durex*np.arange(5))
ax6.text(xtik[0] + args.durex/2., -0.5, '(0, 0)', horizontalalignment='center', verticalalignment='top')
ax6.text(xtik[1] + args.durex/2., -0.5, '(1, 0)', horizontalalignment='center', verticalalignment='top')
ax6.text(xtik[2] + args.durex/2., -0.5, '(0, 1)', horizontalalignment='center', verticalalignment='top')
ax6.text(xtik[3] + args.durex/2., -0.5, '(1, 1)', horizontalalignment='center', verticalalignment='top')
ax6.set_ylabel('ER, hidden 1 [Hz]')
ax6.set_ylim([0, 11])
#ax6.set_title('Hidden 1 before/after learning', fontdict=None, loc='center', fontsize=9)
ax6.set_xticks(xtik)
ax6.set_yticks([0, 10])
ax6.spines['top'].set_visible(False)
ax6.spines['right'].set_visible(False)
ax6.legend(loc='best')
remove_xticklabel(ax6)
#set_axis_color(ax6, custom_colors['bluish_green'])


# hidden 2 before/after
ax7 = fig.add_subplot(339)
ax7.plot(hidden2_rates[before_learning[0]:before_learning[1], 0], hidden2_rates[before_learning[0]:before_learning[1], 2],
         linestyle='dashed', color=custom_colors['blue'], label='before')
ax7.plot(hidden2_rates[before_learning[0]+1:before_learning[1], 0], hidden2_rates[after_learning[0]+1:after_learning[1], 2],
         color=custom_colors['blue'], label='after')
ax7.text(xtik[0] + args.durex/2., -0.5, '(0, 0)', horizontalalignment='center', verticalalignment='top')
ax7.text(xtik[1] + args.durex/2., -0.5, '(1, 0)', horizontalalignment='center', verticalalignment='top')
ax7.text(xtik[2] + args.durex/2., -0.5, '(0, 1)', horizontalalignment='center', verticalalignment='top')
ax7.text(xtik[3] + args.durex/2., -0.5, '(1, 1)', horizontalalignment='center', verticalalignment='top')
ax7.set_ylabel('ER, hidden 2 [Hz]')
ax7.set_ylim([0, 11])
#ax7.set_title('Hidden 2 before/after learning', fontdict=None, loc='center', fontsize=9)
ax7.set_xticks(xtik)
ax7.set_yticks([0, 10])
ax7.spines['top'].set_visible(False)
ax7.spines['right'].set_visible(False)
ax7.legend(loc='best')
remove_xticklabel(ax7)
#set_axis_color(ax7, custom_colors['bluish_green'])


plt.tight_layout()
plt.savefig('../../results/xor/XOR.pdf')
plt.close()
