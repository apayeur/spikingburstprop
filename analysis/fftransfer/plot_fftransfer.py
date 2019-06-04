import sys
sys.path.append('../')
import raster_analysis as ra
import numpy as np
import matplotlib.pyplot as plt
from utils_plotting import custom_colors
from scipy.optimize import curve_fit
import seaborn as sns
plt.style.use('../thesis_mplrc.dms')

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-filesuffix", type=str, help="suffix indicating feedback weights", default="")
parser.add_argument("-resultdir", type=str, help="parent result directory", default="../../results/fftransfer")
parser.add_argument("-datadir", type=str, help="parent data directory", default="../../data/fftransfer")
args = parser.parse_args()


def sigmoid_fit(x, yshift, sigmax, slope, xshift):
    return yshift + sigmax/(1 + np.exp(-(slope*(x - xshift))))


def exp_fit(x, a, slope):
    return a*np.exp(slope*x)


def combined_sig(x, p1, p2):
    A = np.exp(-p2[2]*(p1[3]-p2[3]))
    print('A = {}'.format(A))
    return p2[0] + p2[1]/(1. + A*(p1[1]/(x - p1[0]) - 1.)**(p2[2]/p1[2]))


def e2_simplified(x, p1, p2):
    A = np.exp(-p2[2]*(p1[3]-p2[3]))
    print('A = {}'.format(A))
    return p2[0] + p2[1]/(1. + A*(p1[1]/(x-p1[0]) - 1.))


def e2_oversimplified(x, p1, p2):
    A = np.exp(-p2[2]*(p1[3]-p2[3]))
    print('A = {}'.format(A))
    return p2[0] + (p2[1]/A/p1[1])*(x-p1[0])


def power_law_fit(x, p1, p2):
    return p2[0]*(x/p1[0])**(p2[1]/p1[1])

currents = np.arange(-200, 500, 100)

# BP, ER, BR for population 1, 2, 3
filenames1 = args.datadir + '/fftransfer.0.ras1'
bp1, er1, br1 = ra.bpercurves_from_raster(filenames1, currents)
#filenames1 = args.datadir + '/fftransfer.0.brate1'
#bp1, er1, br1 = ra.bpercurves(filenames1, currents)

filenames2 = args.datadir + '/fftransfer.0.ras2'
bp2, er2, br2 = ra.bpercurves_from_raster(filenames2, currents)
#filenames2 = args.datadir + '/fftransfer.0.brate2'
#bp2, er2, br2 = ra.bpercurves(filenames2, currents)

#filenames3 = args.datadir + '/fftransfer.0.brate3'
#bp3, er3, br3 = ra.bpercurves(filenames3, currents)

#filename_pv = args.datadir + '/fftransfer.0.pv2rate'
#rate_pv =  ra.ficurves(filename_pv, currents)


# Fit event rates to sigmoid
popt1, pcov1 = curve_fit(sigmoid_fit, currents, er1, p0=[0., 15., 0.1, 200.])
popt2, pcov2 = curve_fit(sigmoid_fit, currents, er2, p0=[4., 6., 0.1, 200.])
#popt3, pcov3 = curve_fit(sigmoid_fit, currents, er3, p0=[5., 6., 0.01, 200.])
#popt1, pcov1 = curve_fit(exp_fit, currents, er1, p0=[1., 0.01])
#popt2, pcov2 = curve_fit(exp_fit, currents, er2, p0=[1., 0.01])

popte2ve1, pcove2ve1 = curve_fit(exp_fit, er1, er2, p0=[5., 0.01])
#popte3ve2, pcove3ve2 = curve_fit(exp_fit, er2, er3, p0=[0.1, 0.01])

# Fit PV rates
# popt_pv, pcov_pv = curve_fit(sigmoid_fit, currents, rate_pv, p0=[0., 30., 0.1, 200.])
#  popt_pv_exp, pcov_pv_exp = curve_fit(exp_fit, er1, rate_pv, p0=[5., 0.01])


print('fit parameters for er1: {}'.format(popt1))
print('fit parameters for er2: {}'.format(popt2))
#print('fit parameters for er3: {}'.format(popt3))
print('fit parameters for er2 versus e1: {}'.format(popte2ve1))

# Plotting
panel_letter_pos = [0.1, 1]
plt.figure(figsize=(5, 3/1.5))
plt.subplot(121)
plt.title(r'\textbf{A} $e_{l} = a_l + c_l \sigma[\alpha_l (I_1 - \gamma_l)]$', loc='left')
#plt.text(0, 1, r'\textbf{A}', horizontalalignment='left', verticalalignment='top', transform=plt.gca().transAxes)
plt.semilogy(currents, er1, color=custom_colors['blue'], label='$e_{1}$')
plt.semilogy(currents, er2, '--', color=custom_colors['blue'], label='$e_{2}$')
#plt.plot(currents, er3, ':', color=custom_colors['blue'], label='$e_{3}$')
plt.semilogy(currents, sigmoid_fit(currents, *popt1), color='grey', lw=1)
plt.semilogy(currents, sigmoid_fit(currents, *popt2), '--', color='grey', lw=1)
#plt.plot(currents, sigmoid_fit(currents, *popt3), ':', color='grey', lw=1)
plt.xlabel('$I_1$ [pA]')
plt.ylabel('$e_l$ [Hz]')
sns.despine()
plt.tight_layout()
plt.legend(fontsize=9)
#plt.savefig(args.resultdir + '/ER1.pdf')
#plt.close()

plt.subplot(122)
plt.title(r'\textbf{B} $e_{2} = f(W_2e_1 + n_2)$', loc='left')
#plt.text(*panel_letter_pos, r'\textbf{C}', horizontalalignment='left', verticalalignment='top', transform=plt.gca().transAxes)
plt.semilogy(er1, er2, color=custom_colors['blue'])
plt.semilogy(er1, exp_fit(np.array(er1), *popte2ve1), ':', color='grey', lw=1, label='exp. fit')
#plt.plot(er1, combined_sig(er1 , popt1, popt2), ':', color='grey', lw=1, label='combined sigmoids')
#plt.plot(er1, exp_fit(np.array(er1), *popte2ve1), '--', color='black', lw=1, label='exp. fit')
#plt.plot(er1, e2_oversimplified(er1, popt1, popt2), '--', color='blue', lw=1, label='linear approx.')
plt.xlabel('$e_1$ [Hz]')
plt.ylabel('$e_2$ [Hz]')
sns.despine()
plt.legend(fontsize=9)

'''
plt.subplot(223)
plt.title(r'\textbf{C} $e_{3} = f(W_3e_2 + n_3)$', loc='left')
#plt.text(*panel_letter_pos, r'\textbf{C}', horizontalalignment='left', verticalalignment='top', transform=plt.gca().transAxes)
plt.plot(er2, er3, color=custom_colors['blue'])
plt.plot(er2, e2_simplified(er2, popt2, popt3), '--', color='grey', lw=1, label='combined sigmoids')
plt.plot(er2, exp_fit(np.array(er2), *popte3ve2), '--', color='black', lw=1, label='exp. fit')
#plt.plot(er1, e2_oversimplified(er1, popt1, popt2), '--', color='blue', lw=1, label='linear approx.')
plt.xlabel('$e_2$ [Hz]')
plt.ylabel('$e_3$ [Hz]')
sns.despine()
plt.legend(fontsize=9)

'''
plt.tight_layout()
plt.savefig(args.resultdir + '/FFTransferFunction.pdf')
plt.close()




'''
plt.figure(figsize=(5, 3/1.5))
plt.subplot(121)
plt.title(r'\textbf{A}', loc='left')
plt.plot(currents, rate_pv, color='black', lw=2)
plt.plot(currents, sigmoid_fit(currents, *popt_pv), '--', color='grey', lw=1, label='sigmoid fit')
#plt.plot(currents, exp_fit(currents, *popt_pv_exp), '--', color='grey', lw=1)
plt.xlabel('$I_1$ [pA]')
plt.ylabel('Rate PV [Hz]')
sns.despine()
plt.tight_layout()
plt.legend(fontsize=9)

plt.subplot(122)
plt.title(r'\textbf{B}', loc='left')
plt.plot(er1, rate_pv, color='black', lw=2)
plt.plot(er1, exp_fit(np.array(er1), *popt_pv_exp), '--', color='grey', label='exp. fit')
plt.xlabel('$e_1$ [Hz]')
plt.ylabel('Rate PV [Hz]')
sns.despine()
plt.legend(fontsize=9)
plt.tight_layout()

plt.savefig(args.resultdir + '/FFTransferFunction_PV.pdf')
plt.close()
'''