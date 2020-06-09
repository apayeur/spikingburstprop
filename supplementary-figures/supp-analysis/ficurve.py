import sys
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.gridspec as gridspec
sys.path.append('../../analysis/')
from utils_plotting import custom_colors
plt.style.use('../../analysis/thesis_mplrc.dms')
from raster_analysis import *
import seaborn as sns
from matplotlib import rc
rc('text', usetex=True)


def extract_bp_and_er(fn):
    """
     Compute BP and ER curve based on sequence of current injections
     :param fn: file name
     :return:
     """
    data = np.loadtxt(fn)
    times = data[:, 0]
    bp = 100. * data[:, 1] / data[:, 2]
    br = data[:, 1]
    er = data[:, 2]
    loc_mean_bp = []
    #loc_std_bp = []
    loc_mean_er = []
    loc_mean_br = []
    for i, dendritic_current in enumerate(dendritic_currents):
        indices = np.where(np.logical_and(times > i + 0.5, times < i + 1.))[0]
        loc_mean_bp.append(np.mean(bp[indices]))
        loc_mean_er.append(np.mean(er[indices]))
        loc_mean_br.append(np.mean(br[indices]))
        #loc_std_bp.append(np.std(bp[indices]))
    return loc_mean_bp, loc_mean_er, loc_mean_br


#Set output path
OUTPUTPATH = '../supp-results/normalization'
if not os.path.exists(OUTPUTPATH):
    os.makedirs(OUTPUTPATH)

#Set dendritic currents according to what was used in the simulation (see sim_normalization_of_bpcurves.cpp)
dendritic_currents = np.arange(-100, 350, 50)

#Set input data
files = ['../supp-data/normalization/somacurrent0.000000_wpyrtopyr0.020000.0.brate',
         '../supp-data/normalization/somacurrent100.000000_wpyrtopyr0.020000.0.brate',
         '../supp-data/normalization/somacurrent200.000000_wpyrtopyr0.020000.0.brate']

#Compute BP, ER, BR
mean_bp = []
mean_er = []
mean_br = []
for f in files:
    (l_mean_bp, l_mean_er, l_mean_br) = extract_bp_and_er(f)
    #times, er, br = extract_event_times_and_burst_times_from_raster(f, binsize=5.e-3, tau=16.e-3)
    #l_mean_er, l_mean_br, l_mean_bp = mean_er_and_br(times, er, br, dendritic_currents)
    mean_bp.append(l_mean_bp)
    mean_er.append(l_mean_er)
    mean_br.append(l_mean_br)

#Plotting
f, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.5,2.5))
alphas = [0.3, 0.65, 1]
for i in range(len(mean_bp)):
    if i==0:
        #ax1.plot(dendritic_currents, mean_bp[i], color=custom_colors['red'], alpha=alphas[i], label='BP, $I_\mathrm{soma} =$ ' + files[i][27:30] + ' pA')
        #ax1.plot(dendritic_currents, mean_er[i], color=custom_colors['blue'], alpha=alphas[i], label='ER, $I_\mathrm{soma} =$ ' + files[i][27:30] + ' pA')
        ax1.plot(dendritic_currents, mean_bp[i], color=custom_colors['red'], alpha=alphas[i], label='BP')
        ax1.plot(dendritic_currents, mean_er[i], color=custom_colors['blue'], alpha=alphas[i], label='ER')
    else:
        ax1.plot(dendritic_currents, mean_bp[i], color=custom_colors['red'], alpha=alphas[i])
        ax1.plot(dendritic_currents, mean_er[i], color=custom_colors['blue'], alpha=alphas[i])

ax1.set_xlabel('$I_\mathrm{dend}$ [pA]')
ax1.set_ylabel('BP [\%], ER [Hz]')
ax1.legend(loc='best')

for i in range(len(mean_bp)):
    if i==0:
        #ax2.plot(dendritic_currents, mean_br[i], color=custom_colors['orange'], alpha=alphas[i], label='BR, $I_\mathrm{soma} =$ ' + files[i][27:30] + ' pA')
        #ax2.plot(dendritic_currents, mean_er[i], color=custom_colors['blue'], alpha=alphas[i], label='ER, $I_\mathrm{soma} =$ ' + files[i][27:30] + ' pA')
        ax2.plot(dendritic_currents, mean_br[i], color=custom_colors['orange'], alpha=alphas[i],
                 label='BR')
        ax2.plot(dendritic_currents, mean_er[i], color=custom_colors['blue'], alpha=alphas[i],
                 label='ER')
    else:
        ax2.plot(dendritic_currents, mean_br[i], color=custom_colors['orange'], alpha=alphas[i])
        ax2.plot(dendritic_currents, mean_er[i], color=custom_colors['blue'], alpha=alphas[i])
#ax2.set_xlim(xmax=100)
#ax2.set_ylim(ymax=50)
ax2.set_xlabel('$I_\mathrm{dend}$ [pA]')
ax2.set_ylabel('BR [Hz], ER [Hz]')
ax2.legend(loc='best')
sns.despine()
plt.tight_layout()
#plt.show()
plt.savefig(OUTPUTPATH + '/BPversusIdend.png')
plt.close()