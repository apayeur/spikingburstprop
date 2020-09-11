import sys
sys.path.append('../')
import numpy as np
import argparse
import matplotlib.pyplot as plt
from utils_plotting import set_axis_color
from utils_plotting import custom_colors
from utils_plotting import remove_xticklabel
plt.style.use('../thesis_mplrc.dms')


parser = argparse.ArgumentParser()
parser.add_argument("-datadir", type=str, help="data directory", default=".")
parser.add_argument("-resultdir", type=str, help="result directory", default=".")
parser.add_argument("-connect_type", type=str, help="connection type", default=".")
parser.add_argument("-outfile_name", type=str, help="name of output file", default=".")
parser.add_argument("-durex", type=float, help="duration of example")
parser.add_argument("-alpha", type=float, help="moving average time scale")

args = parser.parse_args()


# --- import data --- #
wsum = np.loadtxt(args.datadir + '/plasticity_rule.0.wsum')              # average of nonzero synaptic weights
brate_post = np.loadtxt(args.datadir + '/plasticity_rule.0.brate')       # ER, BR and BP of postsynaptic population
brate_pre = np.loadtxt(args.datadir + '/plasticity_rule.0.brate_input')  # ER, BR and BP of presynaptic population

# construct average of moving average BP
er_trace = np.loadtxt(args.datadir + '/plasticity_rule0.0.trevent')
br_trace = np.loadtxt(args.datadir + '/plasticity_rule0.0.trburst')
for i in range(1, 50):
    er_tmp = np.loadtxt(args.datadir + '/plasticity_rule'+str(i)+'.0.trevent')
    br_tmp = np.loadtxt(args.datadir + '/plasticity_rule'+str(i)+'.0.trburst')
    er_trace[:, 1] += er_tmp[:, 1]
    br_trace[:, 1] += br_tmp[:, 1]
er_trace[:, 1] /= 50.
br_trace[:, 1] /= 50.


# --- plotting --- #
fig = plt.figure(figsize=(70./25.4, 88./25.4))
tics = list(4*args.alpha + args.durex*np.arange(0, 7))
binSize = brate_post[1, 0] - brate_post[0, 0]

ax1 = fig.add_subplot(311)
ax1.plot(brate_post[:, 0], brate_post[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax1.plot(er_trace[:, 0], np.ones(er_trace[:, 0].shape)*er_trace[int((3*args.alpha + args.durex)/binSize), 1]/args.alpha, ':', color="black", lw=0.5)
ax1.set_ylabel('Postsyn. ER [Hz]', color=custom_colors['blue'])
ax1.tick_params('y', labelcolor=custom_colors['blue'])
ax1.set_ylim([4,9])
ax1.spines['top'].set_visible(False)
ax1.set_xticks([0]+tics)
remove_xticklabel(ax1)
#set_axis_color(ax1, custom_colors['sky_blue'])

ax1_tw = ax1.twinx()
ax1_tw.plot(brate_post[:, 0], 100*brate_post[:, 1]/brate_post[:, 2], color=custom_colors['red'], lw=1)
ax1_tw.plot(er_trace[:, 0], 100*br_trace[:, 1]/er_trace[:, 1], '--', color=custom_colors['red'], lw=1, label=r'$\overline{P}$')
ax1_tw.tick_params('y', labelcolor=custom_colors['red'])
ax1_tw.set_yticks([0, 25, 50])
ax1_tw.set_ylim([0, 60])
ax1_tw.legend(loc='upper right', bbox_to_anchor=(1, 1.1))
ax1_tw.set_ylabel('Postsyn. BP [%]', color=custom_colors['red'])
ax1_tw.spines['top'].set_visible(False)
#set_axis_color(ax1_tw, custom_colors['sky_blue'])

ax2 = fig.add_subplot(312, sharex=ax1)
ax2.plot(wsum[:, 0], 100*(wsum[:, 1] - wsum[0, 1])/wsum[0, 1], color='black', lw=1.5, label='Weight change')
ax2.plot(wsum[:, 0], 100*(brate_post[:, 1]/brate_post[:, 2] - br_trace[:, 1]/er_trace[:, 1]), '-.', lw=1,
         color=custom_colors['red'], label=r'Postsyn. BP$ - \overline{P}$')
ax2.set_ylabel('[%]')
ax2.legend(loc='upper right', bbox_to_anchor=(1, 1.05))
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
remove_xticklabel(ax2)

ax3 = fig.add_subplot(313)
ax3.plot(brate_pre[:, 0], brate_pre[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax3.set_ylim([0,5])
ax3.set_xlabel('Time [s]')
ax3.set_ylabel('Presyn. ER [Hz]', color=custom_colors['blue'])
ax3.tick_params('y', labelcolor=custom_colors['blue'])
ax3.set_xticks([0]+tics)
ax3.spines['top'].set_visible(False)
#set_axis_color(ax3, custom_colors['bluish_green'])

ax3_tw = ax3.twinx()
ax3_tw.set_yticks([0, 25, 50])
ax3_tw.set_ylim([0, 60])
ax3_tw.plot(brate_pre[:, 0], 100*brate_pre[:, 1]/brate_pre[:, 2], color=custom_colors['red'], lw=1)
ax3_tw.set_ylabel('Presyn. BP [%]', color=custom_colors['red'])
ax3_tw.tick_params('y', labelcolor=custom_colors['red'])
ax3_tw.spines['top'].set_visible(False)
#set_axis_color(ax3_tw, custom_colors['bluish_green'])

plt.tight_layout()
plt.savefig(args.outfile_name)
plt.close()
