import sys
sys.path.append('../')
import numpy as np
import argparse
import matplotlib.pyplot as plt
from utils_plotting import set_axis_color
from utils_plotting import custom_colors
from utils_plotting import remove_xticklabel
plt.style.use('../thesis_mplrc.dms')

# Parse arguments from the command line
parser = argparse.ArgumentParser()
parser.add_argument("-datadir", type=str, help="data directory", default=".")
parser.add_argument("-resultdir", type=str, help="result directory", default=".")
parser.add_argument("-connect_type", type=str, help="connection type", default=".")
parser.add_argument("-outfile_name", type=str, help="name of output file", default=".")
parser.add_argument("-durex", type=float, help="duration of example")
parser.add_argument("-alpha", type=float, help="moving average time scale")
parser.add_argument("-numberofrealiz", type=int, help="number of realization of the noise", default=1)

args = parser.parse_args()
number_of_realizations = args.numberofrealiz

# ~~~~~~~~~~~~~~~~~~~~~ #
#       Import data     #
# ~~~~~~~~~~~~~~~~~~~~~ #
# Wsum (average of nonzero synaptic weights
wsum = np.loadtxt('{}/plasticity_rule.0.wsum_seed{}'.format(args.datadir, 1))
for i in range(1, number_of_realizations):
    filename = '{}/plasticity_rule.0.wsum_seed{}'.format(args.datadir, i+1)
    wsum = np.hstack([wsum, np.reshape(np.loadtxt(filename)[:, 1], (-1,1))])  # average of nonzero synaptic weights
wsum_mean = np.mean(wsum[:, 1:], axis=1)
wsum_std = np.std(wsum[:, 1:], axis=1)

# Event and burst rate of the pre and postsynaptic neurons
brate_post = np.loadtxt('{}/plasticity_rule.0.brate_seed{}'.format(args.datadir, 1))
for i in range(1, number_of_realizations):
    filename = '{}/plasticity_rule.0.brate_seed{}'.format(args.datadir, i+1)
    brate_post = np.hstack([brate_post, np.loadtxt(filename)[:, 1:]])  # BR and ER of postsynaptic population
brate_post_meanBR, brate_post_meanER = np.mean(brate_post[:, 1::2], axis=1), np.mean(brate_post[:, 2::2], axis=1)
brate_post_stdBR, brate_post_stdER = np.std(brate_post[:, 1::2], axis=1), np.std(brate_post[:, 2::2], axis=1)
brate_post_meanBP = np.mean(brate_post[:, 1::2]/brate_post[:, 2::2], axis=1)
brate_post_stdBP = np.std(brate_post[:, 1::2]/brate_post[:, 2::2], axis=1)

brate_pre = np.loadtxt('{}/plasticity_rule.0.brate_input_seed{}'.format(args.datadir, 1))
for i in range(1, number_of_realizations):
    filename = '{}/plasticity_rule.0.brate_input_seed{}'.format(args.datadir, i+1)
    brate_pre = np.hstack([brate_pre, np.loadtxt(filename)[:, 1:]])  # BR and ER of presynaptic population
brate_pre_meanBR = np.mean(brate_pre[:, 1::2], axis=1)
brate_pre_meanER = np.mean(brate_pre[:, 2::2], axis=1)
brate_pre_meanBP = np.mean(brate_pre[:, 1::2]/brate_pre[:, 2::2], axis=1)
brate_pre_stdBR = np.std(brate_pre[:, 1::2], axis=1)
brate_pre_stdER = np.std(brate_pre[:, 2::2], axis=1)
brate_pre_stdBP = np.std(brate_pre[:, 1::2]/brate_pre[:, 2::2], axis=1)

# construct average of moving average BP
all_er_traces = []
all_br_traces = []
for r in range(1, number_of_realizations+1):
    er_trace = np.loadtxt(args.datadir + '/plasticity_rule0.0.trevent_seed' + str(r))
    br_trace = np.loadtxt(args.datadir + '/plasticity_rule0.0.trburst_seed' + str(r))
    for i in range(1, 50):
        er_tmp = np.loadtxt(args.datadir + '/plasticity_rule'+str(i)+'.0.trevent_seed'+str(r))
        br_tmp = np.loadtxt(args.datadir + '/plasticity_rule'+str(i)+'.0.trburst_seed'+str(r))
        er_trace[:, 1] += er_tmp[:, 1]
        br_trace[:, 1] += br_tmp[:, 1]
    all_er_traces.append(er_trace[:, 1] / 50.)
    all_br_traces.append(br_trace[:, 1] / 50.)
all_er_traces = np.array(all_er_traces)
all_br_traces = np.array(all_br_traces)
mean_ER_trace = np.mean(all_er_traces, axis=0)
mean_BR_trace = np.mean(all_br_traces, axis=0)
mean_BP_trace = np.mean(all_br_traces/all_er_traces, axis=0)
std_BP_trace = np.std(all_br_traces/all_er_traces, axis=0)

# --- plotting --- #
fig = plt.figure(figsize=(70./25.4, 88./25.4))
tics = list(4*args.alpha + args.durex*np.arange(0, 7))
binSize = brate_post[1, 0] - brate_post[0, 0]

ax1 = fig.add_subplot(311)
ax1.plot(brate_post[:, 0], brate_post_meanER, color=custom_colors['blue'], lw=1, label='ER')
ax1.fill_between(brate_post[:, 0],
                 brate_post_meanER - 2*brate_post_stdER,
                 brate_post_meanER + 2*brate_post_stdER,
                 color=custom_colors['blue'], alpha=0.5, lw=0)
ax1.plot(brate_post[:, 0], np.ones_like(brate_post_meanER)*np.mean(brate_post_meanER[int(args.alpha/binSize):int((3*args.alpha + args.durex)/binSize)]), ':', color="black", lw=0.5)
ax1.set_ylabel('Postsyn. ER [Hz]', color=custom_colors['blue'])
ax1.set_ylim([0, 9])
ax1.set_yticks([0, 8])
ax1.spines['top'].set_visible(False)
ax1.set_xticks([0]+tics)
remove_xticklabel(ax1)

ax1_tw = ax1.twinx()
ax1_tw.plot(brate_post[:, 0], 100*brate_post_meanBP, color=custom_colors['red'], lw=1)
ax1_tw.fill_between(brate_post[:, 0],
                    100*brate_post_meanBP - 2*100*brate_post_stdBP,
                    100*brate_post_meanBP + 2*100*brate_post_stdBP,
                    color=custom_colors['red'], alpha=0.5, lw=0)
ax1_tw.plot(er_trace[:, 0], 100*mean_BR_trace/mean_ER_trace, '--', color=custom_colors['red'], lw=1, label=r'$\overline{P}$')
ax1_tw.fill_between(er_trace[:, 0],
                    100*mean_BP_trace - 2*100*std_BP_trace,
                    100*mean_BP_trace + 2*100*std_BP_trace,
                    color=custom_colors['red'], alpha=0.5, lw=0)
ax1_tw.set_yticks([0, 50])
ax1_tw.set_ylim([0, 60])
ax1_tw.legend(loc='upper right', bbox_to_anchor=(1, 1.1))
ax1_tw.set_ylabel('Postsyn. BP [%]', color=custom_colors['red'])
ax1_tw.spines['top'].set_visible(False)

ax2 = fig.add_subplot(312, sharex=ax1)
ax2.plot(wsum[:, 0], 100*(wsum_mean - wsum_mean[0])/wsum[0, 1], color='black', lw=1.5, label='Weight change')
ax2.fill_between(wsum[:, 0], 100*(wsum_mean - wsum_mean[0])/wsum[0, 1] - 100*2*wsum_std/wsum[0, 1],
                 100*(wsum_mean - wsum_mean[0])/wsum[0, 1] + 100*2*wsum_std/wsum[0, 1],
                 color='black', alpha=0.5, lw=0)
ax2.plot(wsum[:, 0], 100*(brate_post_meanBP - mean_BP_trace), '-.', lw=1,
         color=custom_colors['red'], label=r'Postsyn. BP$ - \overline{P}$')
ax2.fill_between(wsum[:, 0],
                 100*(brate_post_meanBP - mean_BP_trace) - 2*100*np.sqrt(brate_post_stdBP**2 + std_BP_trace**2),
                 100*(brate_post_meanBP - mean_BP_trace) + 2*100*np.sqrt(brate_post_stdBP**2 + std_BP_trace**2),
                 color=custom_colors['red'], alpha=0.5, lw=0)
ax2.plot(wsum[:, 0], np.zeros_like(wsum[:, 0]), ':', color="black", lw=0.5)
ax2.set_ylabel('[%]')
ax2.legend(loc='upper right', bbox_to_anchor=(1, 1.05))
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
remove_xticklabel(ax2)

ax3 = fig.add_subplot(313)
ax3.plot(brate_pre[:, 0], brate_pre_meanER, color=custom_colors['blue'], lw=1, label='ER')
ax3.fill_between(brate_pre[:, 0],
                 brate_pre_meanER - 2*brate_pre_stdER,
                 brate_pre_meanER + 2*brate_pre_stdER,
                 color=custom_colors['blue'], alpha=0.5, lw=0)
ax3.set_ylim([0, 6])
ax3.set_yticks([0, 5])
ax3.set_xlabel('Time [s]')
ax3.set_ylabel('Presyn. ER [Hz]', color=custom_colors['blue'])
ax3.set_xticks([0]+tics)
ax3.spines['top'].set_visible(False)

ax3_tw = ax3.twinx()
ax3_tw.set_yticks([0, 50])
ax3_tw.set_ylim([0, 60])
ax3_tw.plot(brate_pre[:, 0], 100*brate_pre_meanBP, color=custom_colors['red'], lw=1)
ax3_tw.fill_between(brate_pre[:, 0],
                    100*brate_pre_meanBP - 2*100*brate_pre_stdBP,
                    100*brate_pre_meanBP + 2*100*brate_pre_stdBP,
                    color=custom_colors['red'], alpha=0.5, lw=0)
ax3_tw.set_ylabel('Presyn. BP [%]', color=custom_colors['red'])
ax3_tw.spines['top'].set_visible(False)

plt.tight_layout()
plt.savefig(args.outfile_name)
plt.close()
