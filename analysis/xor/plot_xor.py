import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import sys
sys.path.append('../')
plt.style.use('../thesis_mplrc.dms')
from utils_plotting import remove_xticklabel
from utils_plotting import custom_colors

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-durex", type=float, help="duration of example")
parser.add_argument("-alpha", type=float, help="moving average time scale")
parser.add_argument("-numex", type=int, help="number of training examples")
args = parser.parse_args()


# import data
data_file_prefix = '../../data/xor/xor.0.brate_'
datadir = '../../data/xor/'

output_rates = np.loadtxt(data_file_prefix + 'output')[-int(4*args.durex/1.):, :]
hidden1_rates = np.loadtxt(data_file_prefix + 'hidden1')[-int(4*args.durex/1.):, :]
hidden2_rates = np.loadtxt(data_file_prefix + 'hidden2')[-int(4*args.durex/1.):, :]
input1_rates = np.loadtxt(data_file_prefix + 'input1')[-int(4*args.durex/1.):, :]
input2_rates = np.loadtxt(data_file_prefix + 'input2')[-int(4*args.durex/1.):, :]
er_trace = np.loadtxt(datadir + 'xor.0.trevent')[-int(4*args.durex/1.):, :]
br_trace = np.loadtxt(datadir + 'xor.0.trburst')[-int(4*args.durex/1.):, :]

data_file_prefix = '../../data/xor/xor.0.wsum_'
wsum_in1_to_hid1 = np.loadtxt(data_file_prefix + 'in1_to_hid1')
wsum_in1_to_hid2 = np.loadtxt(data_file_prefix + 'in1_to_hid2')
wsum_in2_to_hid1 = np.loadtxt(data_file_prefix + 'in2_to_hid1')
wsum_in2_to_hid2 = np.loadtxt(data_file_prefix + 'in2_to_hid2')
wsum_hid1_to_out = np.loadtxt(data_file_prefix + 'hid1_to_out')
wsum_hid2_to_out = np.loadtxt(data_file_prefix + 'hid2_to_out')

# position of xticks according with protocol
#tics = list(3*args.alpha + args.durex*np.arange(0, args.numex+1))

# plotting rates and burst probabilities
fig = plt.figure(1, figsize=(6, 6/1.5))
gs1 = gridspec.GridSpec(3, 2)

#  output
ax1 = fig.add_subplot(gs1[0, 1])
#ax1.set_aspect(1/1.5, adjustable='box', anchor='C')
ax1.plot(output_rates[:, 0], output_rates[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax1.plot(output_rates[:, 0], output_rates[:, 1], color=custom_colors['orange'], lw=1, label='BR')
ax1.set_ylabel('Output rates [Hz]')
ax1.set_ylim([0,15])
#ax1.set_xlim([3006, 3026])
#ax1.set_xticks([0]+tics)
ax1.spines['top'].set_visible(False)
ax1.legend(loc='best')

ax1_tw = ax1.twinx()
ax1_tw.plot(output_rates[:, 0], 100*output_rates[:, 1]/output_rates[:, 2], color=custom_colors['red'], lw=1)
ax1_tw.plot(er_trace[:, 0], 100*br_trace[:, 1]/er_trace[:, 1], '--', color=custom_colors['red'], lw=1)
ax1_tw.tick_params('y', labelcolor=custom_colors['red'])
ax1_tw.set_yticks([0, 30, 60, 90])
#ax1_tw.set_xlim([3006, 3026])
ax1_tw.set_ylabel('Postsyn. BP [\%]', color=custom_colors['red'])
remove_xticklabel(ax1)


#  hidden layer 1
ax2 = fig.add_subplot(gs1[1, 0])
ax2.plot(hidden1_rates[:, 0], hidden1_rates[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax2.plot(hidden1_rates[:, 0], hidden1_rates[:, 1], color=custom_colors['orange'], lw=1, label='BR')
ax2.set_ylabel('Hidden1 rates [Hz]')
#ax2.set_xticks([0]+tics)
#ax2.set_xlim([3006, 3026])
ax2.set_ylim([0,15])
ax2.spines['top'].set_visible(False)

ax2_tw = ax2.twinx()
ax2_tw.plot(hidden1_rates[:, 0], 100*hidden1_rates[:, 1]/hidden1_rates[:, 2], color=custom_colors['red'], lw=1)
ax2_tw.tick_params('y', labelcolor=custom_colors['red'])
ax2_tw.set_yticks([0, 30, 60, 90])
#ax2_tw.set_xlim([3006, 3026])
ax2_tw.set_ylabel('Postsyn. BP [\%]', color=custom_colors['red'])
remove_xticklabel(ax2)

#  hidden layer 2
ax3 = fig.add_subplot(gs1[1, 1])
ax3.plot(hidden2_rates[:, 0], hidden2_rates[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax3.plot(hidden2_rates[:, 0], hidden2_rates[:, 1], color=custom_colors['orange'], lw=1, label='BR')
ax3.set_ylabel('Hidden2 rates [Hz]')
#ax3.set_xlim([3006, 3026])
ax3.set_ylim([0,15])
#ax3.set_xticks([0]+tics)
ax3.spines['top'].set_visible(False)

ax3_tw = ax3.twinx()
ax3_tw.plot(hidden2_rates[:, 0], 100*hidden2_rates[:, 1]/hidden2_rates[:, 2], color=custom_colors['red'], lw=1)
ax3_tw.tick_params('y', labelcolor=custom_colors['red'])
ax3_tw.set_yticks([0, 30, 60, 90])
#ax3_tw.set_xlim([3006, 3026])
ax3_tw.set_ylabel('Postsyn. BP [\%]', color=custom_colors['red'])
remove_xticklabel(ax3)

#  input layer 1
ax4 = fig.add_subplot(gs1[2, 0])
ax4.plot(input1_rates[:, 0], input1_rates[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax4.set_xlabel("Time [s]")
ax4.set_ylabel('ER, input1 [Hz]')
#ax4.set_xlim([3006, 3026])
ax4.set_ylim([0,15])
#ax4.set_xticks([0]+tics)
#ax4.spines['top'].set_visible(False)

#  input layer 2
ax5 = fig.add_subplot(gs1[2, 1])
ax5.plot(input2_rates[:, 0], input2_rates[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax5.set_xlabel("Time [s]")
ax5.set_ylabel('ER, input2 [Hz]')
#ax5.set_xlim([3006, 3026])
ax5.set_ylim([0,15])
#ax5.set_xticks([0]+tics)
#ax5.spines['top'].set_visible(False)

plt.tight_layout()
plt.savefig('../../results/xor/Activity.pdf')
plt.close()


# plotting sum of weights
fig = plt.figure(2, figsize=(5, 5/1.5))
ax1 = fig.add_subplot(221)
ax1.plot(wsum_hid1_to_out[:,0], wsum_hid1_to_out[:,1])
ax1.set_ylabel('$\Sigma$ weights')
ax1.set_title('Hidden 1 to output')
#ax1.set_xticks([0]+tics)
remove_xticklabel(ax1)

ax2 = fig.add_subplot(222)
ax2.plot(wsum_hid2_to_out[:,0], wsum_hid2_to_out[:,1])
ax2.set_title('Hidden 2 to output')
#ax2.set_xticks([0]+tics)
remove_xticklabel(ax2)

ax3 = fig.add_subplot(223)
ax3.plot(wsum_in1_to_hid1[:,0], wsum_in1_to_hid1[:,1], label='to hidden 1')
ax3.plot(wsum_in1_to_hid2[:,0], wsum_in1_to_hid2[:,1], label='to hidden 2')
ax3.set_title('Input 1 to hidden layer')
ax3.set_xlabel('Time [s]')
ax3.set_ylabel('$\Sigma$ weights')
#ax3.set_xticks([0]+tics)
ax3.legend()

ax4 = fig.add_subplot(224)
ax4.plot(wsum_in2_to_hid1[:,0], wsum_in2_to_hid1[:,1], label='to hidden 1')
ax4.plot(wsum_in2_to_hid2[:,0], wsum_in2_to_hid2[:,1], label='to hidden 2')
ax4.set_title('Input 2 to hidden layer')
ax4.set_xlabel('Time [s]')
#ax4.set_xticks([0]+tics)
ax4.legend()

plt.tight_layout()
plt.savefig('../../results/xor/Wsum.pdf')
plt.close()