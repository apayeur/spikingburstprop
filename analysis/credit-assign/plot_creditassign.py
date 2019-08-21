import matplotlib.pyplot as plt
import sys
sys.path.append('../')
import raster_analysis as ra
import voltage_analysis as va
import numpy as np
import argparse
import utils_plotting as up
plt.style.use('../thesis_mplrc.dms')

from build_current import example_credit_assign

parser = argparse.ArgumentParser()
parser.add_argument("-filesuffix", type=str, help="suffix indicating FF weights", default="")
parser.add_argument("-resultdir", type=str, help="parent result directory", default="../../results/credit-assign")

args = parser.parse_args()

# load data
brate1 = np.loadtxt('../../data/credit-assign/credit_assign.0.brate1_seed'+str(1))
brate2 = np.loadtxt('../../data/credit-assign/credit_assign.0.brate2_seed'+str(1))

er_trace = np.loadtxt('../../data/credit-assign/credit_assign0.0.trevent')
br_trace = np.loadtxt('../../data/credit-assign/credit_assign0.0.trburst')
for i in range(1,50):
    er_tmp = np.loadtxt('../../data/credit-assign/credit_assign'+str(i)+'.0.trevent')
    br_tmp = np.loadtxt('../../data/credit-assign/credit_assign'+str(i)+'.0.trburst')
    er_trace[:,1] += er_tmp[:,1]
    br_trace[:,1] += br_tmp[:,1]
er_trace[:,1] /= 50.
br_trace[:,1] /= 50.

wsum = np.loadtxt('../../data/credit-assign/credit_assign.0.wsum')


# import currents
current_soma = np.loadtxt("../../data/credit-assign/current_soma.txt")
current_dend = np.loadtxt("../../data/credit-assign/current_dend.txt")
binSize = current_soma[1,0] - current_soma[0,0]


# rescale BR2
BP1 = 100*brate1[:,1]/brate1[:,2]
min_BP1 = np.min(BP1[int(1/binSize):])
max_BP1 = np.max(BP1[int(1/binSize):])
BR2 = brate2[:,1]
min_BR2 = np.min(BR2[int(1/binSize):])
max_BR2 = np.max(BR2[int(1/binSize):])
rescaled_BR2 = min_BP1 + (BR2 - min_BR2)*(max_BP1-min_BP1)/(max_BR2-min_BR2)

# plotting
plt.figure(figsize=(4, 4/1.6))
tics = [0, 10, 25, 40, 55]
ax = plt.subplot(221)
plt.plot(brate1[:,0], BP1, color=up.custom_colors['red'])
plt.plot(brate1[:,0], rescaled_BR2, '--', color=up.custom_colors['orange'], lw=1, label='BR2 (scaled)')
plt.xlim(xmin=0)
plt.ylim([0., 100.])
plt.xticks(tics)
plt.ylabel(r'BP [\%]')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend(loc='best')

plt.subplot(222, sharex=ax)
plt.plot(brate2[:,0], 100*brate2[:,1]/brate2[:,2], lw=1, color=up.custom_colors['red'])
plt.plot(er_trace[:,0], 100*br_trace[:,1]/er_trace[:,1], '--', color=up.custom_colors['red'], lw=1, label='Reference BP')
plt.xlim(xmin=0)
plt.ylim([0., 100.])
plt.ylabel(r'BP [\%]')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.legend(loc='best')

plt.subplot(223, sharex=ax)
plt.plot(brate1[:,0], brate1[:,2], color=up.custom_colors['blue'])
plt.xlim(xmin=0)
plt.ylim([0., 15])
plt.xlabel('Time [s]')
plt.ylabel('ER [Hz]')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.legend(loc='best')

plt.subplot(224, sharex=ax)
plt.plot(brate2[:,0], brate2[:,2], color=up.custom_colors['blue'])
plt.xlim(xmin=0)
plt.ylim([0., 15])
plt.xlabel('Time [s]')
plt.ylabel('ER [Hz]')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.legend(loc='best')

plt.tight_layout()
plt.savefig(args.resultdir+'/CreditAssign' + args.filesuffix + '.pdf')
plt.close()


plt.figure(figsize=(2.5, 2.5/1.6))
plt.plot(wsum[:, 0], wsum[:, 1], color='black', lw=1.5)
plt.xticks(tics)
plt.xlabel('Time [s]')
plt.ylabel('Average weight')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(args.resultdir+'/Weight' + args.filesuffix + '.pdf')
plt.close()