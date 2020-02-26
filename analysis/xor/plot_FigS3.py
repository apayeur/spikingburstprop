import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import sys
sys.path.append('../')
plt.style.use('../thesis_mplrc.dms')
plt.rcParams["axes.labelsize"] = 9
from utils_plotting import remove_xticklabel
from utils_plotting import custom_colors
import seaborn as sns


# data folders
datadir = '../../data/xor/'
folder_names = ['very-good-xor-corrected',
                'very-good-xor-corrected-finite-size',
                'very-good-xor-corrected-no-hidden-plasticity',
                'very-good-xor-corrected-random-examples']

# import costs
cost = dict()
for folder_name in folder_names:
    cost[folder_name] = np.loadtxt(datadir + folder_name + '/cost_epoch.dat')

# import output ERs
outputER = dict()
for folder_name in folder_names:
    output_rate = np.loadtxt(datadir + folder_name + '/xor.0.brate_output')
    # select "after learning" data (20 seconds/example * 4 examples)
    binSize = output_rate[1, 0] - output_rate[0, 0]
    after_learning = [output_rate.shape[0]-int(4*20./binSize), output_rate.shape[0]]
    outputER[folder_name] = output_rate[output_rate.shape[0]-int(4*20./binSize):, 2]

# compute threshold
e_min = 2.
e_max = 10.
e_th = e_min + 0.5*(e_max - e_min)

# plotting
fig = plt.figure(figsize=(5, 3))
gs = gridspec.GridSpec(2,4)
ax = dict()
ax_costs = fig.add_subplot(gs[:2, :2])
ax_normal = fig.add_subplot(gs[0, 2])
ax_no_hidden_plasticity = fig.add_subplot(gs[0, 3])
ax_finite_size = fig.add_subplot(gs[1, 2])
ax_random_examples = fig.add_subplot(gs[1, 3])

colors_sns = sns.color_palette('colorblind')
linet = ['-', '--', ':', '-.']

# plot costs
ax_costs.plot(cost[folder_names[0]][:,0]+1, cost[folder_names[0]][:,1],
              linet[0], color=colors_sns[0], lw=1, label='Fig.4')
ax_costs.plot(cost[folder_names[1]][:,0]+1, cost[folder_names[1]][:,1],
              linet[1], color=colors_sns[1], lw=1, label='400 neurons')
ax_costs.plot(cost[folder_names[2]][:,0]+1, cost[folder_names[2]][:,1],
              linet[2], color=colors_sns[2], lw=1, label='No hidden plasticity')
ax_costs.plot(cost[folder_names[3]][:,0]+1, cost[folder_names[3]][:,1],
              linet[3], color=colors_sns[3], lw=1, label='Randomized examples')
ax_costs.annotate('A', (-0.2, 1), xycoords='axes fraction', weight='bold')
ax_costs.set_xlabel('Epoch')
ax_costs.set_ylabel('Cost')
ax_costs.set_xticks([1, 500])
ax_costs.set_yticks([2., 2.6, 3.2])
sns.despine(ax=ax_costs, trim=True)
ax_costs.legend()

# plot output ERs
ax_normal.plot(outputER[folder_names[0]], linet[0], color=colors_sns[0])
ax_normal.plot(e_th*np.ones(outputER[folder_names[0]].shape), '--', color='grey', lw=1)
ax_normal.set_ylabel('Output ER [Hz]')
ax_normal.set_xticks([], [])
ax_normal.annotate(r'\textbf{B}', (-0.5, 1), xycoords='axes fraction')
sns.despine(ax=ax_normal, bottom=True, trim=True)

ax_finite_size.plot(outputER[folder_names[1]], linet[1], color=colors_sns[1])
ax_finite_size.plot(e_th*np.ones(outputER[folder_names[1]].shape), '--', color='grey', lw=1)
ax_finite_size.set_ylabel('Output ER [Hz]')
ax_finite_size.set_xticks([], [])
ax_finite_size.annotate(r'\textbf{D}', (-0.5,1), xycoords='axes fraction')
sns.despine(ax=ax_finite_size, bottom=True, trim=True)

ax_no_hidden_plasticity.plot(outputER[folder_names[2]], linet[2], color=colors_sns[2])
ax_no_hidden_plasticity.plot(e_th*np.ones(outputER[folder_names[2]].shape), '--', color='grey', lw=1)
ax_no_hidden_plasticity.set_xticks([], [])
ax_no_hidden_plasticity.annotate(r'\textbf{C}', (-0.4,1), xycoords='axes fraction')
sns.despine(ax=ax_no_hidden_plasticity, bottom=True, trim=True)

ax_random_examples.plot(outputER[folder_names[3]], linet[3], color=colors_sns[3])
ax_random_examples.plot(e_th*np.ones(outputER[folder_names[3]].shape), '--', color='grey', lw=1)
ax_random_examples.set_xticks([], [])
ax_random_examples.annotate(r'\textbf{E}', (-0.4,1), xycoords='axes fraction')
sns.despine(ax=ax_random_examples, bottom=True, trim=True)

plt.tight_layout()
plt.savefig('../../results/xor/FigureS3.png')
plt.close()