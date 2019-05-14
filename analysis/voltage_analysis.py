import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import sys
import seaborn as sns
#sys.path.append('./')
from utils_plotting import custom_colors
plt.style.use('../thesis_mplrc.dms')


def output_histogram(filename, output_file, comp_type='soma'):
    V = np.loadtxt(filename)
    V = 1000*V[int(1./1e-4):, 1]
    plt.figure(figsize=(2.*1.6, 2.))
    plt.hist(V[V < -40], 500, normed=True)
    if comp_type == 'soma':
        threshold = -45
    elif comp_type == 'dend':
        threshold = -60
    meanV = np.mean(V[V < threshold])
    stdV  = np.std(V[V < threshold])

    plt.text(0.9, 0.8, 'mean = {:2.0f}'.format(meanV), horizontalalignment='right', verticalalignment = 'center', transform = plt.gca().transAxes)
    plt.text(0.9, 0.7, 'std  = {:3.1f}'.format(stdV), horizontalalignment='right', verticalalignment = 'center', transform = plt.gca().transAxes)
    plt.xlim(xmax=-45)
    plt.xlabel('Voltage [mV]')
    plt.ylabel('PDF')
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


def population_voltage(filenames, threshold=-40e-3):
    """
    Compute the population-average voltage, using V < threshold.
    :param filenames: names of files containing voltage data
    :param threshold:
    :return: time and voltage arrays
    """
    for i, filename in enumerate(filenames):
        data = np.loadtxt(filename)
        V_tmp = data[:, 1]
        mask_tmp = np.ones(V_tmp.shape)
        indices = np.where(V_tmp > threshold)[0]
        V_tmp[indices] = 0.
        mask_tmp[indices] = 0.
        if i == 0:
            t = data[:, 0]
            V = V_tmp
            mask = mask_tmp
        else:
            V += V_tmp
            mask += mask_tmp

    V = V/mask

    return t, V

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