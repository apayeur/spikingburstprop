import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import sys
#sys.path.append('./')
from utils_plotting import custom_colors
plt.style.use('../analysis/thesis_mplrc.dms')


def rates_from_raster(filename_raster, binsize=20.e-3, tau=16.e-3):
    """
    :param filename_raster: raster file (1st column = times, 2nd  column = neuron number)
    :param binsize: Size of the discretization bins (in seconds)
    :param tau: burst detection threshold (in seconds)
    :return:
    """
    raster = np.loadtxt(filename_raster)
    nb_of_bins = int(np.max(raster[:, 0]) / binsize)
    nb_of_neurons = int(max(raster[:, 1]))+1
    event_rate = np.zeros(nb_of_bins+1)
    burst_rate = np.zeros(nb_of_bins+1)

    for n in range(nb_of_neurons):
        spiketimes = raster[np.where(raster[:, 1] == n)[0], 0]
        if len(spiketimes) > 0:
            events = [spiketimes[0]]
            event_rate[int(spiketimes[0] / binsize)] += 1.
            bursts = []

            burst_state = -1
            for s in range(1, len(spiketimes)):
                if spiketimes[s] - spiketimes[s-1] > tau:
                    events.append(spiketimes[s])
                    burst_state=-1.
                    event_rate[int(spiketimes[s]/binsize)] += 1.
                else:
                    if burst_state < 0:
                        bursts.append(spiketimes[s-1]) # s-1 to correct for the fact that bursts are detected on the second spike
                        burst_state=1
                        burst_rate[int(spiketimes[s-1] / binsize)] += 1.

    event_rate = event_rate/(nb_of_neurons * binsize)
    burst_rate = burst_rate/(nb_of_neurons * binsize)
    times = binsize*np.arange(0, len(event_rate))
    times += binsize/2.
    return times, event_rate, burst_rate


def rates_from_brates(filename):
    data = np.loadtxt(filename)
    BR = data[:,1]
    ER = data[:,2]
    return data[:,0],ER, BR


def average_rates_over_realizations(filenames, binsize=20.e-3, tau=16.e-3):
    """
    Computes the sample average event and burst rates
    :param filenames:   List of raster files
    :param binsize:     Size of the discretization bins (in seconds)
    :param tau:         Burst detection threshold (in seconds)
    :return:            Time-dependent averages event and burst rates, over different seeds for the random noise
    """
    mean_ER = np.array([])
    mean_BR = np.array([])
    for i, filename in enumerate(filenames):
        times, ER, BR = rates_from_raster(filename, binsize=binsize, tau=tau)
        #times, ER, BR = rates_from_brates(filename)
        if i==0:
            mean_ER = ER
            mean_BR = BR
        else:
            mean_ER += ER
            mean_BR += BR
    mean_ER /= len(filenames)
    mean_BR /= len(filenames)
    return times, mean_ER, mean_BR


def plot_rates(filenames, binsize=20.e-3, tau=16.e-3):
    """
    Display average rates and burst probability.
    :param filenames:   List of raster files
    :param binsize:     Size of the discretization bins (in seconds)
    :param tau:         Burst detection threshold (in seconds)
    """
    t, ER, BR = average_rates_over_realizations(filenames, binsize=binsize, tau=tau)
    plt.plot(t, ER, color=custom_colors['blue'], label='ER')
    plt.plot(t, BR, color=custom_colors['orange'], label='BR')
    plt.plot(t, 100*BR/ER, color=custom_colors['red'], label='BP')
    plt.xlim([1,3])
    plt.xlabel('Time [ms]')
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_rates_with_inputs(filenames, inputs, binsize=20.e-3, tau=16.e-3):
    """
     Display average rates, burst probability together with the dendritic and somatic inputs
     :param filenames:   List of raster files
     :param inputs:      Dictionary with element 'dendrite' and 'soma' and 'times'
     :param binsize:     Size of the discretization bins (in seconds)current traces in pA
     :param tau:         Burst detection threshold (in seconds)
     """
    t, ER, BR = average_rates_over_realizations(filenames, binsize=binsize, tau=tau)
    BP = 100*BR/ER
    plt.figure(figsize=(5, 5))
    plt.subplot(311)
    plt.plot(t, BP, color=custom_colors['red'], label='BP')
    rescaled_input = np.min(BP) + (inputs['dendrite'] - np.min(inputs['dendrite'])) *\
                                   (np.max(BP)-np.min(BP))/(np.max(inputs['dendrite'])-np.min(inputs['dendrite']))
    plt.plot(inputs['times'], rescaled_input, '--', color=custom_colors['red'])
    plt.xlim([inputs['times'][0], inputs['times'][-1]])
    plt.ylabel(r'BP, $I_\mathrm{d}$ (scaled)')

    plt.subplot(312)
    plt.plot(t, BR, color=custom_colors['orange'], label='BR')
    plt.xlim([inputs['times'][0], inputs['times'][-1]])
    plt.ylabel('BR [Hz]')

    plt.subplot(313)
    plt.plot(t, ER, color=custom_colors['blue'], label='ER')
    rescaled_input = np.min(ER) + (inputs['soma'] - np.min(inputs['soma'])) *\
                                   (np.max(ER)-np.min(ER))/(np.max(inputs['soma'])-np.min(inputs['soma']))
    plt.plot(inputs['times'], rescaled_input, '--', color=custom_colors['blue'])
    plt.xlim([inputs['times'][0], inputs['times'][-1]])
    plt.xlabel('Time [s]')
    plt.ylabel('ER [Hz], $I_\mathrm{s}$ (scaled)')
    plt.tight_layout()
    plt.savefig('MultiplexingWithBurstPoisson.pdf')
    plt.close()
    #plt.show()


def add_step(bg, begin, length, amplitude):
    '''
    :param bg: background np.array
    :param begin: beginning of step (index)
    :param length: duration (# of bins)
    :param amplitude: amplitude of the step relative to bg
    :return:
    '''
    if begin + length < len(bg):
        bg[begin:begin+length+1] += amplitude


def alternating_square_current(times, mini, maxi, begins, duration):
    """
    Construct an alternating square current with minimum current *min* and maximum current *max*.
    :param times:   Array of times
    :param min:     Minimum current amplitude
    :param max:     Maximum current amplitude
    :param begins:  Array of times at which the current jumps up
    :param duration: Duration of square (in second)
    :return:
    """
    current = mini*np.ones(times.shape)
    dt = times[1] - times[0]
    for begin in begins:
        add_step(current, int(begin/dt), int(duration/dt), maxi - mini)

    return current


def mean_er_and_br(times, er, br, dendritic_currents):
    loc_mean_er = []
    loc_mean_br = []
    loc_mean_bp = []

    for i, dendritic_current in enumerate(dendritic_currents):
        indices = np.where(np.logical_and(times > i + 0.5, times < i + 1.))[0]
        loc_mean_er.append(np.mean(er[indices]))
        loc_mean_br.append(np.mean(br[indices]))
        loc_mean_bp.append(np.mean(100*br[indices]/er[indices]))
    return np.array(loc_mean_er), np.array(loc_mean_br), np.array(loc_mean_bp)



'''
if __name__ == "__main__":
    dendritic_currents = np.arange(-100, 250, 25.)
    filename = 'ficurve.0.ras_somacurrent100.000000'
    t, ER, BR = er_and_br_from_raster(filename, binsize=10e-3, tau=16e-3)
    mean_ER, mean_BR = mean_er_and_br(t, ER, BR, dendritic_currents)
    plt.plot(dendritic_currents, mean_ER, color=custom_colors['blue'])
    plt.plot(dendritic_currents, 100*mean_BR/mean_ER, color=custom_colors['red'])
    #plt.plot(t, ER, label='ER')
    #plt.plot(t, BR/ER*100, label='BP')
    #plt.plot(t, BR, label='BR')
    plt.legend()
    plt.tight_layout()
    plt.show()
'''