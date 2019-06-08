import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import sys
import seaborn as sns
import raster_analysis as ra
#sys.path.append('./')
from utils_plotting import custom_colors


def example_credit_assign(params, binsize=5.e-3):
    '''
    Build the currents injected into the somas of population 1
    and the dendrites of population 2.
    :return: times, somatic current [A], dendritic current [A]
    '''
    baseline_soma = params['baseline_soma']
    baseline_dend = params['baseline_dend']
    example_duration = int(params['example_duration']/binsize)  # duration of input applied to soma (# of bins)
    relaxation_period = int(params['relaxation_period']/binsize)
    nb_of_ex = params['number_of_examples']
    amplitudes_soma = params['amplitudes_soma']  # list of amplitudes
    amplitudes_dend = params['amplitudes_dend']  # list of amplitudes
    slope_duration = int(params['slope_duration']/binsize)

    assert nb_of_ex == len(amplitudes_soma)
    assert nb_of_ex == len(amplitudes_dend)

    # Construct time array (initial relaxation = 2*relaxation_period
    total_duration = 2 * relaxation_period + nb_of_ex * (example_duration + relaxation_period)
    t = binsize * np.arange(0, total_duration)

    # Construct currents
    current_soma = baseline_soma*np.ones(t.shape)
    for i in range(nb_of_ex):
        ra.add_soft_step(current_soma, 2*relaxation_period + i*(example_duration+relaxation_period),
                         example_duration, amplitudes_soma[i], slope_duration)

    current_dend = baseline_dend*np.ones(t.shape)
    for i in range(nb_of_ex):
        ra.add_step(current_dend, 2*relaxation_period + example_duration//2 + i*(example_duration+relaxation_period),
                         example_duration//2, amplitudes_dend[i])


    np.savetxt("../data/credit-assign/current_soma.txt", (np.array([t, current_soma])).T, fmt=['%e', '%e'])
    np.savetxt("../data/credit-assign/current_dend.txt", (np.array([t, current_dend])).T, fmt=['%e', '%e'])

    return t, current_soma, current_dend


if __name__ == '__main__':
    # Testing function example_credit_assign
    params = {'baseline_soma': 100.e-12,
              'baseline_dend': 50.e-12,
              'example_duration': 500.e-3,
              'relaxation_period': 500.e-3,
              'number_of_examples': 3,
              'amplitudes_soma': [100.e-12, 100.e-12, -50.e-12],
              'amplitudes_dend': [50.e-12, -50.e-12, 50.e-12],
              'slope_duration': 0.1*500.e-3}

    t, current_soma, current_dend = example_credit_assign(params, binsize=20.e-3)

    plt.figure(figsize=(5,5/1.6))
    plt.plot(t, current_soma/1.e-12, '--', color=custom_colors['blue'], label='soma current')
    plt.plot(t, current_dend/1.e-12, '--', color=custom_colors['red'], label='dend. current')
    plt.xlabel('Time [s]')
    plt.ylabel('Current [pA]')
    plt.legend()
    plt.tight_layout()
    plt.show()






