import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../')

from utils_plotting import remove_xticklabel
from utils_plotting import custom_colors
plt.style.use('../thesis_mplrc.dms')

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-datadir", type=str, help="data directory", default=".")
parser.add_argument("-resultdir", type=str, help="result directory", default=".")
parser.add_argument("-connect_type", type=str, help="connection type", default=".")
parser.add_argument("-outfile_name", type=str, help="name of output file", default=".")

args = parser.parse_args()


# import data
wsum = np.loadtxt(args.datadir + '/plasticity_rule.0.wsum')
brate_post = np.loadtxt(args.datadir + '/plasticity_rule.0.brate')
brate_pre = np.loadtxt(args.datadir + '/plasticity_rule.0.brate_input')


# plotting
fig = plt.figure(figsize=(4, 7./1.6))

ax1 = fig.add_subplot(311)
ax1.plot(brate_post[:, 0], brate_post[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax1.plot(brate_post[:, 0], brate_post[:, 1], color=custom_colors['orange'], lw=1, label='BR')
ax1.set_ylabel('Postsyn. rates [Hz]')
ax1.spines['top'].set_visible(False)
ax1.legend(loc='best')
#ax1.set_xticks([0, 200, 400, 600, 800, 1000, 1200])
remove_xticklabel(ax1)

ax1_tw = ax1.twinx()
ax1_tw.plot(brate_post[:, 0], 100*brate_post[:, 1]/brate_post[:, 2], color=custom_colors['red'], lw=1)
ax1_tw.tick_params('y', labelcolor=custom_colors['red'])
ax1_tw.set_yticks([0, 30, 60, 90])
ax1_tw.set_ylabel('Postsyn. BP [\%]', color=custom_colors['red'])
ax1_tw.spines['top'].set_visible(False)

ax2 = fig.add_subplot(312, sharex=ax1)
ax2.plot(wsum[:, 0], wsum[:, 1], color='black', lw=1.5)
ax2.set_ylabel('$\sum$ somatic weights')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
#ax2.set_xticks([0, 200, 400, 600, 800, 1000, 1200])
remove_xticklabel(ax2)

ax3 = fig.add_subplot(313)
ax3.plot(brate_pre[:, 0], brate_pre[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax3.plot(brate_pre[:, 0], brate_pre[:, 1], color=custom_colors['orange'], lw=1, label='BR')
ax3.set_ylim(ymax=5)
ax3.set_xlabel('Time [s]')
ax3.set_ylabel('Presyn. rates [Hz]')
#ax3.set_xticks([0, 200, 400, 600, 800, 1000, 1200])
ax3.spines['top'].set_visible(False)

ax3_tw = ax3.twinx()
ax3_tw.plot(brate_pre[:, 0], 100*brate_pre[:, 1]/brate_pre[:, 2], color=custom_colors['red'], lw=1)
ax3_tw.set_ylabel('Presyn. BP [\%]', color=custom_colors['red'])
ax3_tw.tick_params('y', labelcolor=custom_colors['red'])
ax3_tw.spines['top'].set_visible(False)

plt.tight_layout()
plt.savefig(args.outfile_name)
plt.close()
