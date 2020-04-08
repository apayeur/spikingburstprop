import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import sys
sys.path.append('../')
plt.style.use('../thesis_mplrc.dms')
from utils_plotting import remove_xticklabel
from utils_plotting import custom_colors
from utils_plotting import set_axis_color
from raster_analysis import add_step
import seaborn as sns

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-durex", type=float, help="duration of example")
parser.add_argument("-alpha", type=float, help="moving average time scale")
parser.add_argument("-numex", type=int, help="number of training examples")
args = parser.parse_args()


def mean_output(segment):
    # select data segments
    r_out = output_rates[segment[0]:segment[1]+1, 2]
    durex_in_bins = int(args.durex/binSize)
    discard = 1 #  discard data points near transitions from one example to the other
    r00 = np.mean(r_out[discard:durex_in_bins-discard])
    r10 = np.mean(r_out[durex_in_bins+discard:2*durex_in_bins-discard])
    r01 = np.mean(r_out[2*durex_in_bins+discard:3*durex_in_bins-discard])
    r11 = np.mean(r_out[3*durex_in_bins+discard:-discard])
    return r00, r10, r01, r11


def std_output(segment):
    # select data segments
    r_out = output_rates[segment[0]:segment[1]+1, 2]
    durex_in_bins = int(args.durex/binSize)
    discard = 1
    r00 = np.std(r_out[discard:durex_in_bins - discard])
    r10 = np.std(r_out[durex_in_bins + discard:2 * durex_in_bins - discard])
    r01 = np.std(r_out[2 * durex_in_bins + discard:3 * durex_in_bins - discard])
    r11 = np.std(r_out[3 * durex_in_bins + discard:-discard])
    return r00, r10, r01, r11


# import data
data_file_prefix = '../../data/xor/very-good-xor-corrected/xor.0.brate_'
datadir = '../../data/xor/very-good-xor-corrected/'

output_rates = np.loadtxt(data_file_prefix + 'output')
hidden1_rates = np.loadtxt(data_file_prefix + 'hidden1')
hidden2_rates = np.loadtxt(data_file_prefix + 'hidden2')
input1_rates = np.loadtxt(data_file_prefix + 'input1')
input2_rates = np.loadtxt(data_file_prefix + 'input2')
cost_epoch = np.loadtxt(datadir + 'cost_epoch.dat')

binSize = output_rates[1,0] - output_rates[0,0]
after_learning = [output_rates.shape[0]-int(4*args.durex/binSize), output_rates.shape[0]+1]
before_learning = [int(3*args.alpha/binSize)+1, int(3*args.alpha/binSize)+int(4*args.durex/binSize)+1]
during_learning = [int(3*args.alpha/binSize)+int(args.numex//2*args.durex/binSize)+1,
                   int(3*args.alpha/binSize)+int(args.numex//2*args.durex/binSize)+int(4*args.durex/binSize)+1]

er_trace = np.loadtxt(datadir + 'xor0.0.trevent')
br_trace = np.loadtxt(datadir + 'xor0.0.trburst')
for i in range(1,50):
    er_tmp = np.loadtxt(datadir + 'xor'+str(i)+'.0.trevent')
    br_tmp = np.loadtxt(datadir + 'xor'+str(i)+'.0.trburst')
    er_trace[:,1] += er_tmp[:,1]
    br_trace[:,1] += br_tmp[:,1]
er_trace[:,1] /= 50.
br_trace[:,1] /= 50.
er_trace = er_trace
br_trace = br_trace
er_trace_ex = np.loadtxt(datadir + 'xor0.0.trevent')

lim_seg = after_learning
r_out = output_rates[lim_seg[0]:lim_seg[1], :]
r_hid1 = hidden1_rates[lim_seg[0]:lim_seg[1], :]
r_hid2 = hidden2_rates[lim_seg[0]:lim_seg[1], :]
r_inp1 = input1_rates[lim_seg[0]:lim_seg[1], :]
r_inp2 = input2_rates[lim_seg[0]:lim_seg[1], :]
er_est = er_trace[lim_seg[0]:lim_seg[1], :]
br_est = br_trace[lim_seg[0]:lim_seg[1], :]

xtik = list(lim_seg[0] + args.durex*np.arange(5)-1)

e_min = 2.
e_max = 10.
e_th = e_min + 0.5*(e_max - e_min)

# compute mean +/- SD
m = mean_output(before_learning)
s = std_output(before_learning)
print('Before learning:')
print('(0,0): {:1.2} +/- {:1.1}'.format(m[0], 2*s[0]))
print('(1,0): {:1.2} +/- {:1.1}'.format(m[1], 2*s[1]))
print('(0,1): {:1.2} +/- {:1.1}'.format(m[2], 2*s[2]))
print('(1,1): {:1.2} +/- {:1.1}'.format(m[3], 2*s[3]))

m = mean_output(after_learning)
s = std_output(after_learning)
print('After learning:')
print('(0,0): {:1.2} +/- {:1.1}'.format(m[0], 2*s[0]))
print('(1,0): {:1.2} +/- {:1.1}'.format(m[1], 2*s[1]))
print('(0,1): {:1.2} +/- {:1.1}'.format(m[2], 2*s[2]))
print('(1,1): {:1.2} +/- {:1.1}'.format(m[3], 2*s[3]))


fig = plt.figure(figsize=(183/25.4, 4.5/7.5*183/25.4))
# first subplot is empty

# inputs
ax1 = fig.add_subplot(337)
ax1.plot(r_inp1[1:, 0], r_inp1[1:, 2], color=custom_colors['blue'], label='input 1')
ax1.plot(r_inp2[1:, 0], r_inp2[1:, 2], ':', lw=2, color=custom_colors['light_blue'], label='input 2')
ax1.text(xtik[0] + args.durex/2., -0.5, '(0, 0)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax1.text(xtik[1] + args.durex/2., -0.5, '(1, 0)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax1.text(xtik[2] + args.durex/2., -0.5, '(0, 1)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax1.text(xtik[3] + args.durex/2., -0.5, '(1, 1)', horizontalalignment='center', verticalalignment='top', fontsize=8)
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
ax2.plot(output_rates[before_learning[0]:before_learning[1], 0], output_rates[before_learning[0]:before_learning[1], 2],
         linestyle='dashed', color=custom_colors['blue'], label='before')
ax2.plot(output_rates[before_learning[0]:before_learning[1], 0], output_rates[after_learning[0]:after_learning[1], 2],
         color=custom_colors['blue'], label='after')
ax2.plot(output_rates[before_learning[0]:before_learning[1], 0],
         e_th*np.ones(np.shape(output_rates[before_learning[0]:before_learning[1], 0])), 'k:', lw=0.5)
ax2.text(before_learning[0] + args.durex/2., -0.5, '(0, 0)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax2.text(before_learning[0] + 3*args.durex/2., -0.5, '(1, 0)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax2.text(before_learning[0] + 5*args.durex/2., -0.5, '(0, 1)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax2.text(before_learning[0] + 7*args.durex/2., -0.5, '(1, 1)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax2.set_ylabel('ER, output [Hz]')
ax2.set_ylim([0, 11])
#ax2.set_title('Output before/after learning', fontdict=None, loc='center', fontsize=9)
ax2.set_xticks(list(before_learning[0] + args.durex*np.arange(5)-1))
ax2.set_yticks([0, 10])
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.legend(loc='best')
remove_xticklabel(ax2)
set_axis_color(ax2, custom_colors['sky_blue'])


# error generation at output layer
ax3 = fig.add_subplot(333)
r_out = output_rates[during_learning[0]:during_learning[1], :]
teacher = np.zeros(np.shape(r_out[:,2]))
er_est = er_trace[during_learning[0]:during_learning[1], 1]/args.alpha
add_step(teacher, int(0.9*args.durex/binSize), int(0.1*args.durex/binSize)-1, -1./(e_max - er_est[int(0.9*args.durex/binSize)]))
add_step(teacher, int(1.9*args.durex/binSize), int(0.1*args.durex/binSize)-1, 1./(er_est[int(1.9*args.durex/binSize)] - e_min))
add_step(teacher, int(2.9*args.durex/binSize), int(0.1*args.durex/binSize)-1, 1./(er_est[int(2.9*args.durex/binSize)] - e_min))
add_step(teacher, int(3.9*args.durex/binSize), int(0.1*args.durex/binSize)-1, -1./(e_max - er_est[int(3.9*args.durex/binSize)]))

ax3.plot(r_out[:, 0], 100*r_out[:, 1]/r_out[:, 2], color=custom_colors['red'])
ax3.text(during_learning[0] + args.durex/2., -4.5, '(0, 0)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax3.text(during_learning[0] + 3*args.durex/2., -4.5, '(1, 0)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax3.text(during_learning[0] + 5*args.durex/2., -4.5, '(0, 1)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax3.text(during_learning[0] + 7*args.durex/2., -4.5, '(1, 1)', horizontalalignment='center', verticalalignment='top', fontsize=8)
#ax3.set_title('Error generation at output layer', fontdict=None, loc='center', fontsize=9)
remove_xticklabel(ax3)
set_axis_color(ax3, custom_colors['sky_blue'])
ax3.set_yticks([0, 30, 60, 90])
ax3.set_ylabel('BP, output [%]', color=custom_colors['red'])
ax3.spines['top'].set_visible(False)

ax3_tw = ax3.twinx()
ax3_tw.plot(r_out[:, 0], teacher, '--', lw=1, color=custom_colors['pink'])
ax3_tw.set_xticks(list(during_learning[0] + args.durex*np.arange(5)-1))
set_axis_color(ax3_tw, custom_colors['sky_blue'])
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
ax4.plot(r_out[:, 0], rescaled_BR, '--', color=custom_colors['orange'], label='output BR (scaled)')
ax4.set_xticks(list(during_learning[0] + args.durex*np.arange(5)-1))
ax4.set_yticks([0, 50])
ax4.set_ylim([0, 52])
ax4.set_ylabel('BP, hidden 1 [%]')
ax4.text(during_learning[0] + args.durex/2., -2.5, '(0, 0)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax4.text(during_learning[0] + 3*args.durex/2., -2.5, '(1, 0)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax4.text(during_learning[0] + 5*args.durex/2., -2.5, '(0, 1)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax4.text(during_learning[0] + 7*args.durex/2., -2.5, '(1, 1)', horizontalalignment='center', verticalalignment='top', fontsize=8)
#ax4.set_title('Error propagation to hidden 1', fontdict=None, loc='center', fontsize=9)
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.legend(loc='best')
remove_xticklabel(ax4)
set_axis_color(ax4, custom_colors['bluish_green'])


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
ax5.plot(r_out[:, 0], rescaled_BR, '--', color=custom_colors['orange'], label='output BR (scaled)')

ax5.set_xticks(list(during_learning[0] + args.durex*np.arange(5)-1))
ax5.set_yticks([0, 50])
ax5.set_ylim([0, 52])
ax5.set_ylabel('BP, hidden 2 [%]')
ax5.text(during_learning[0] + args.durex/2., -2.5, '(0, 0)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax5.text(during_learning[0] + 3*args.durex/2., -2.5, '(1, 0)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax5.text(during_learning[0] + 5*args.durex/2., -2.5, '(0, 1)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax5.text(during_learning[0] + 7*args.durex/2., -2.5, '(1, 1)', horizontalalignment='center', verticalalignment='top', fontsize=8)
#ax5.set_title('Error propagation to hidden 2', fontdict=None, loc='center', fontsize=9)
ax5.spines['top'].set_visible(False)
ax5.spines['right'].set_visible(False)
ax5.legend(loc='lower left')
remove_xticklabel(ax5)
set_axis_color(ax5, custom_colors['bluish_green'])

# hidden 1 before/after
ax6 = fig.add_subplot(338)
ax6.plot(hidden1_rates[before_learning[0]:before_learning[1], 0], hidden1_rates[before_learning[0]:before_learning[1], 2],
         linestyle='dashed', color=custom_colors['blue'], label='before')
ax6.plot(hidden1_rates[before_learning[0]+1:before_learning[1], 0], hidden1_rates[after_learning[0]+1:after_learning[1], 2],
         color=custom_colors['blue'], label='after')
ax6.text(before_learning[0] + args.durex/2., -0.5, '(0, 0)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax6.text(before_learning[0] + 3*args.durex/2., -0.5, '(1, 0)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax6.text(before_learning[0] + 5*args.durex/2., -0.5, '(0, 1)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax6.text(before_learning[0] + 7*args.durex/2., -0.5, '(1, 1)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax6.set_ylabel('ER, hidden 1 [Hz]')
ax6.set_ylim([0, 11])
#ax6.set_title('Hidden 1 before/after learning', fontdict=None, loc='center', fontsize=9)
ax6.set_xticks(list(before_learning[0] + args.durex*np.arange(5)-1))
ax6.set_yticks([0, 10])
ax6.spines['top'].set_visible(False)
ax6.spines['right'].set_visible(False)
ax6.legend(loc='best')
remove_xticklabel(ax6)
set_axis_color(ax6, custom_colors['bluish_green'])


# hidden 2 before/after
ax7 = fig.add_subplot(339)
ax7.plot(hidden2_rates[before_learning[0]:before_learning[1], 0], hidden2_rates[before_learning[0]:before_learning[1], 2],
         linestyle='dashed', color=custom_colors['blue'], label='before')
ax7.plot(hidden2_rates[before_learning[0]+1:before_learning[1], 0], hidden2_rates[after_learning[0]+1:after_learning[1], 2],
         color=custom_colors['blue'], label='after')
ax7.text(before_learning[0] + args.durex/2., -0.5, '(0, 0)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax7.text(before_learning[0] + 3*args.durex/2., -0.5, '(1, 0)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax7.text(before_learning[0] + 5*args.durex/2., -0.5, '(0, 1)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax7.text(before_learning[0] + 7*args.durex/2., -0.5, '(1, 1)', horizontalalignment='center', verticalalignment='top', fontsize=8)
ax7.set_ylabel('ER, hidden 2 [Hz]')
ax7.set_ylim([0, 11])
#ax7.set_title('Hidden 2 before/after learning', fontdict=None, loc='center', fontsize=9)
ax7.set_xticks(list(before_learning[0] + args.durex*np.arange(5)-1))
ax7.set_yticks([0, 10])
ax7.spines['top'].set_visible(False)
ax7.spines['right'].set_visible(False)
ax7.legend(loc='best')
remove_xticklabel(ax7)
set_axis_color(ax7, custom_colors['bluish_green'])


plt.tight_layout()
plt.savefig('../../results/xor/XOR_NEW.pdf')
plt.close()
