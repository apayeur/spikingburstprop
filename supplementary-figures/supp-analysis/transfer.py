import matplotlib.pyplot as plt
from matplotlib.cbook import get_sample_data
import numpy as np
import matplotlib.gridspec as gridspec
import sys
import os
sys.path.append('../../analysis/')
from utils_plotting import custom_colors
plt.style.use('../../analysis/thesis_mplrc.dms')
from matplotlib import rc
rc('text', usetex=True)

def extract_bp_curve(fn):
    """
    Compute BP curve based on sequence of current injections
    :param fn: file name
    :return:
    """
    data = np.loadtxt(fn)
    bp = 100. * data[:, 1] / data[:, 2]
    times = data[:, 0]
    loc_mean_bp = []
    loc_std_bp = []
    for i, dendritic_current in enumerate(dendritic_currents):
        indices = np.where(np.logical_and(times > i + 0.5, times < i + 1.))[0]
        loc_mean_bp.append(np.mean(bp[indices]))
        loc_std_bp.append(np.std(bp[indices]))
    return loc_mean_bp, loc_std_bp

#Set output path
OUTPUTPATH = '../supp-results/dendritic-modulation'
if not os.path.exists(OUTPUTPATH):
    os.makedirs(OUTPUTPATH)

#Set types of modulation
modulations = {'closeloop' : ['openloop', 'closeloop_wsomtopyr0_0005',
                         'closeloop_wsomtopyr0_001'],
                'disinhibition_through_vip' : ['closeloop_wsomtopyr0_001',
                         'vipdisinhibition_50',
                         'vipdisinhibition_100'],
               'releaseprob': ['releaseprob0_015', 'closeloop_wsomtopyr0_001', 'releaseprob0_025'],
               'burstiness' : ['closeloop_wsomtopyr0_001', 'burstiness1300', 'burstiness1400']}

#Set dendritic currents according to what was used in the simulation (see sim_dendritic_transfer_func_modulation.cpp)
dendritic_currents = np.linspace(-50, 500, int((500 + 50) / 50) + 1)

#Compute BP curves
mean_bp = dict()
std_bp = dict()
for mod in modulations.keys():
    mean_bp[mod] = []
    std_bp[mod] = []
    for submod in modulations[mod]:
        filename = '../supp-data/dendritic-modulation/' + submod + '.0.brate'
        (l_mean_bp, l_std_bp) = extract_bp_curve(filename)
        mean_bp[mod].append(l_mean_bp)
        std_bp[mod].append(l_std_bp)
    mean_bp[mod] = np.array(mean_bp[mod])
    std_bp[mod] = np.array(std_bp[mod])


# Plotting all panels
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(7,7./1.5))
colors = [custom_colors['red'], custom_colors['bluish_green'], custom_colors['pink']]

alphas = np.linspace(0.3, 1, 3)

axes[0, 0].plot(dendritic_currents, mean_bp['closeloop'][0, :], '--', color='black', lw=2, label='openloop')
for i in range(1, mean_bp['closeloop'].shape[0]):
    axes[0, 0].plot(dendritic_currents, mean_bp['closeloop'][i, :], '-', color=custom_colors['red'], lw=2, alpha=alphas[i-1])
axes[0,0].set_xlim([-50,400])
axes[0,0].set_ylabel('Burst probability [\%]')
axes[0,0].spines['top'].set_visible(False)
axes[0,0].spines['right'].set_visible(False)
axes[0,0].set_title(r'\textbf{A} SOM-to-PC synaptic strength', loc='left')
axes[0,0].legend(loc='lower right')
#newax = fig.add_axes([0.1, 0.75, 0.17, 0.17], anchor='SW', zorder=2)
#newax.imshow(im)
#newax.axis('off')

for i in range(0, mean_bp['disinhibition_through_vip'].shape[0]):
    axes[0, 1].plot(dendritic_currents, mean_bp['disinhibition_through_vip'][i, :], '-', color=custom_colors['bluish_green'], lw=2, alpha=alphas[i])
axes[0,1].spines['top'].set_visible(False)
axes[0,1].spines['right'].set_visible(False)
axes[0,1].set_title(r'\textbf{B} VIP-mediated disinhibition', loc='left')
#newax = fig.add_axes([0.55, 0.75, 0.17, 0.17], anchor='SW', zorder=2)
#newax.imshow(im2)
#newax.axis('off')

alphas = np.linspace(0.2, 1, 4)
for i in range(0, mean_bp['releaseprob'].shape[0]):
    axes[1, 0].plot(dendritic_currents, mean_bp['releaseprob'][i, :], '-', color=custom_colors['blue'], lw=2, alpha=alphas[i])
#newax = fig.add_axes([0.1, 0.3, 0.17, 0.17], anchor='SW', zorder=2)
#newax.imshow(im3)
#newax.axis('off')
axes[1,0].spines['top'].set_visible(False)
axes[1,0].spines['right'].set_visible(False)
axes[1,0].set_title(r'\textbf{C} Release probability PC-to-SOM', loc='left')
axes[1,0].set_xlabel('Dendritic input [pA]')
axes[1,0].set_ylabel('Burst probability [\%]')


alphas = np.linspace(0.3, 1, 3)
for i in range(0, mean_bp['burstiness'].shape[0]):
    axes[1, 1].plot(dendritic_currents, mean_bp['burstiness'][i, :], '-', color=custom_colors['orange'], lw=2, alpha=alphas[i])
#newax = fig.add_axes([0.55, 0.3, 0.17, 0.17], anchor='SW', zorder=2)
#newax.imshow(im4)
#newax.axis('off')
axes[1,1].spines['top'].set_visible(False)
axes[1,1].spines['right'].set_visible(False)
axes[1,1].set_title(r'\textbf{D} Dendritic excitability', loc='left')
axes[1,1].set_xlabel('Dendritic input [pA]')

plt.tight_layout()
plt.savefig(OUTPUTPATH+'/DendriticTransfer.pdf')
plt.close()
#plt.show()
