import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import sys
import seaborn as sns
#sys.path.append('./')
from utils_plotting import custom_colors
plt.style.use('thesis_mplrc.dms')


def output_histogram(filename):
    







if __name__ == '__main__':
    #  import voltage traces
    Vm1 = np.loadtxt('../data/propagation/propagation.0.memsoma1')
    Vm2 = np.loadtxt('../data/propagation/propagation.0.memsoma2')

    Vd1 = np.loadtxt('../data/propagation/propagation.0.Vd')
    Vd2 = np.loadtxt('../data/propagation/propagation.0.Vd2')

    VmPV  = np.loadtxt('../data/propagation/propagation.0.mempv')
    VmSOM = np.loadtxt('../data/propagation/propagation.0.memsom')

    #  plotting
    fig = plt.figure(figsize=(5,5))

    fig.add_subplot(321)
    plt.hist(1000*Vm1[int(1./1e-4):, 1], 500)
    plt.xlim(xmax=-45)
    plt.xlabel('Voltage [mV]')
    plt.ylabel('PDF (soma 1)')

    fig.add_subplot(322)
    plt.hist(1000*Vm2[int(1./1e-4):, 1], 500)
    plt.xlim(xmax=-45)
    plt.xlabel('Voltage [mV]')
    plt.ylabel('PDF (soma 2)')

    fig.add_subplot(323)
    plt.hist(1000*Vd1[int(1./1e-4):, 1], 500)
    plt.xlabel('Voltage [mV]')
    plt.ylabel('PDF (dend. 1)')

    fig.add_subplot(324)
    plt.hist(1000*Vd2[int(1./1e-4):, 1], 500)
    plt.xlabel('Voltage [mV]')
    plt.ylabel('PDF (dend. 2)')

    fig.add_subplot(325)
    plt.hist(1000*VmPV[int(1./1e-4):, 1], 500)
    plt.xlim(xmax=-30)
    plt.xlabel('Voltage [mV]')
    plt.ylabel('PDF (PV)')

    fig.add_subplot(326)
    plt.hist(1000*VmSOM[int(1./1e-4):, 1], 500)
    plt.xlim(xmax=-30)
    plt.xlabel('Voltage [mV]')
    plt.ylabel('PDF (SOM)')

    plt.tight_layout()
    plt.show()