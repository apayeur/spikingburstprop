import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import sys
sys.path.append('../')
plt.style.use('../thesis_mplrc.dms')
plt.rcParams["axes.labelsize"] = 9
import seaborn as sns


# data folders
datadir = '../../data/xor-test/'
folder_names = ['sym-ff/multiple-reals',
                'alpha2_durex4',
                'alpha0.5_durex8']

example_durations = {folder_names[0]: 8.,
                     folder_names[1]: 4.,
                     folder_names[2]: 8.}

# import costs
cost = dict()
for folder_name in folder_names:
    cost[folder_name] = np.loadtxt(datadir + folder_name + '/cost_epoch.dat')

# import output ERs
outputER = dict()
for folder_name in folder_names:
    output_rate = np.loadtxt(datadir + folder_name + '/xor-test.0.brate_output_seed1')
    binSize = output_rate[1, 0] - output_rate[0, 0]
    after_learning = [output_rate.shape[0]-int(4*example_durations[folder_name]/binSize), output_rate.shape[0]]
    outputER[folder_name] = output_rate[output_rate.shape[0]-int(4*example_durations[folder_name]/binSize):, 2]

# threshold
e_min = 2.
e_max = 10.
e_th = e_min + 0.5*(e_max - e_min)

# plotting
fig, (ax_cost, ax_output) = plt.subplots(nrows=1, ncols=2, figsize=(5, 2.5))
colors_sns = sns.color_palette('colorblind')
linet = ['-', '--', ':', '-.']

# plot costs
ax_cost.plot(cost[folder_names[0]][:,0]+1, cost[folder_names[0]][:,1],
              linet[0], color=colors_sns[0], lw=1, label=r'$\tau_\mathrm{avg}=2.0 \mathrm{~s}, T=8 \mathrm{~s}$')
ax_cost.plot(cost[folder_names[1]][:,0]+1, cost[folder_names[1]][:,1],
              linet[1], color=colors_sns[1], lw=1, label=r'$\tau_\mathrm{avg}=2.0 \mathrm{~s}, T=4 \mathrm{~s}$')
ax_cost.plot(cost[folder_names[2]][:,0]+1, cost[folder_names[2]][:,1],
              linet[2], color=colors_sns[2], lw=1, label=r'$\tau_\mathrm{avg}=0.5 \mathrm{~s}, T=8 \mathrm{~s}$')

ax_cost.annotate('a', (-0.2, 1), xycoords='axes fraction', weight='bold')
ax_cost.set_xlabel('Epoch')
ax_cost.set_ylabel('Cost')
ax_cost.set_xticks([1, 500])
ax_cost.set_yticks([1., 2., 3.])
sns.despine(ax=ax_cost, trim=True)
ax_cost.legend()

# plot output ERs
for i in range(3):
    ax_output.plot(outputER[folder_names[i]], linet[i], color=colors_sns[i])
ax_output.plot(e_th*np.ones(outputER[folder_names[0]].shape), '--', color='grey', lw=1)
ax_output.set_ylabel('Output ER [Hz]')
ax_output.set_xticks([], [])
ax_output.annotate('b', (-0.5, 1), xycoords='axes fraction')
sns.despine(ax=ax_output, bottom=True, trim=True)


plt.tight_layout()
plt.savefig('../../results/xor-test/FigureSupp_TimeScale.pdf')
plt.close()