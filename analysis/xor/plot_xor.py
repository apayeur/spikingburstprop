import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import sys
sys.path.append('../')
plt.style.use('../thesis_mplrc.dms')
from utils_plotting import custom_colors

# import data
data_file_prefix = '../../data/xor/xor.0.brate_'

output_rates = np.loadtxt(data_file_prefix + 'output')
hidden1_rates = np.loadtxt(data_file_prefix + 'hidden1')
hidden2_rates = np.loadtxt(data_file_prefix + 'hidden2')
input1_rates = np.loadtxt(data_file_prefix + 'input1')
input2_rates = np.loadtxt(data_file_prefix + 'input2')

#'../../results/credit-assign/PV'

# plotting rates and burst probabilities
fig = plt.figure(1, figsize=(6, 6))
gs1 = gridspec.GridSpec(3, 2)

#  output
ax1 = fig.add_subplot(gs1[0, :])
ax1.set_aspect('equal', adjustable='box', anchor='C')
ax1.plot(output_rates[:, 0], output_rates[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax1.plot(output_rates[:, 0], output_rates[:, 1], color=custom_colors['orange'], lw=1, label='BR')
ax1.set_xlabel("Time [s]")
ax1.set_ylabel('Output rates [Hz]')
ax1.spines['top'].set_visible(False)
ax1.legend(loc='best')

ax1_tw = ax1.twinx()
ax1_tw.plot(output_rates[:, 0], 100*output_rates[:, 1]/output_rates[:, 2], color=custom_colors['red'], lw=1)
ax1_tw.tick_params('y', labelcolor=custom_colors['red'])
ax1_tw.set_yticks([0, 30, 60, 90])
ax1_tw.set_ylabel('Postsyn. BP [\%]', color=custom_colors['red'])

#  hidden layer 1
ax2 = fig.add_subplot(gs1[1, 0])
ax2.plot(hidden1_rates[:, 0], hidden1_rates[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax2.plot(hidden1_rates[:, 0], hidden1_rates[:, 1], color=custom_colors['orange'], lw=1, label='BR')
ax2.set_xlabel("Time [s]")
ax2.set_ylabel('Hidden1 rates [Hz]')
ax2.spines['top'].set_visible(False)

ax2_tw = ax2.twinx()
ax2_tw.plot(hidden1_rates[:, 0], 100*hidden1_rates[:, 1]/hidden1_rates[:, 2], color=custom_colors['red'], lw=1)
ax2_tw.tick_params('y', labelcolor=custom_colors['red'])
ax2_tw.set_yticks([0, 30, 60, 90])
ax2_tw.set_ylabel('Postsyn. BP [\%]', color=custom_colors['red'])

#  hidden layer 2
ax3 = fig.add_subplot(gs1[1, 1])
ax3.plot(hidden2_rates[:, 0], hidden2_rates[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax3.plot(hidden2_rates[:, 0], hidden2_rates[:, 1], color=custom_colors['orange'], lw=1, label='BR')
ax3.set_xlabel("Time [s]")
ax3.set_ylabel('Hidden2 rates [Hz]')
ax3.spines['top'].set_visible(False)

ax3_tw = ax3.twinx()
ax3_tw.plot(hidden2_rates[:, 0], 100*hidden2_rates[:, 1]/hidden2_rates[:, 2], color=custom_colors['red'], lw=1)
ax3_tw.tick_params('y', labelcolor=custom_colors['red'])
ax3_tw.set_yticks([0, 30, 60, 90])
ax3_tw.set_ylabel('Postsyn. BP [\%]', color=custom_colors['red'])

#  input layer 1
ax4 = fig.add_subplot(gs1[2, 0])
ax4.plot(input1_rates[:, 0], input1_rates[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax4.plot(input1_rates[:, 0], input1_rates[:, 1], color=custom_colors['orange'], lw=1, label='BR')
ax4.set_xlabel("Time [s]")
ax4.set_ylabel('Input1 rates [Hz]')
ax4.spines['top'].set_visible(False)

ax4_tw = ax4.twinx()
ax4_tw.plot(input1_rates[:, 0], 100*input1_rates[:, 1]/input1_rates[:, 2], color=custom_colors['red'], lw=1)
ax4_tw.tick_params('y', labelcolor=custom_colors['red'])
ax4_tw.set_yticks([0, 30, 60, 90])
ax4_tw.set_ylabel('Postsyn. BP [\%]', color=custom_colors['red'])

#  input layer 2
ax5 = fig.add_subplot(gs1[2, 1])
ax5.plot(input2_rates[:, 0], input2_rates[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax5.plot(input2_rates[:, 0], input2_rates[:, 1], color=custom_colors['orange'], lw=1, label='BR')
ax5.set_xlabel("Time [s]")
ax5.set_ylabel('Input2 rates [Hz]')
ax5.spines['top'].set_visible(False)

ax5_tw = ax5.twinx()
ax5_tw.plot(input2_rates[:, 0], 100*input2_rates[:, 1]/input2_rates[:, 2], color=custom_colors['red'], lw=1)
ax5_tw.tick_params('y', labelcolor=custom_colors['red'])
ax5_tw.set_yticks([0, 30, 60, 90])
ax5_tw.set_ylabel('Postsyn. BP [\%]', color=custom_colors['red'])

plt.tight_layout()
plt.savefig('../../results/xor/Activity.pdf')
plt.close()