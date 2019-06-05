import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import sys
import seaborn as sns
#sys.path.append('./')
from utils_plotting import custom_colors
plt.style.use('../thesis_mplrc.dms')


#################################################################
#           FUNCTIONS RELATED TO RATE COMPUTATIONS              #
#################################################################
def BRER_from_raster(filename_raster, binsize=20.e-3, tau=16.e-3):
    """
    Computes the event rate (ER) and burst rate (BR) from a raster file.
    :param filename_raster: Raster file name (1st column = times, 2nd  column = neuron number)
    :param binsize:         Size of the discretization bins (in seconds)
    :param tau:             Burst detection threshold (in seconds)
    :returns:
        - times - centers of bins
        - event_rate - event rates in Hz
        - burst_rate - burst rates in Hz
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
    times += binsize/2.  # center of bin
    return times, event_rate, burst_rate


def rates_from_raster(filename_raster, binsize=20.e-3):
    """
    Computes rate (PSTH) from a raster file.
    :param filename_raster: Raster file name (1st column = times, 2nd  column = neuron number)
    :param binsize:         Size of the discretization bins (in seconds)
    :returns:
        - times - centers of bins
        - rate  - PSTH in Hz
    """
    raster = np.loadtxt(filename_raster)
    nb_of_bins = int(np.max(raster[:, 0]) / binsize)
    nb_of_neurons = int(max(raster[:, 1]))+1
    rate = np.zeros(nb_of_bins+1)

    for spiketime in raster[:, 0]:
        rate[int(spiketime/binsize)] += 1.

    rate = rate/(nb_of_neurons * binsize)
    times = binsize*np.arange(0, len(rate))
    times += binsize/2.  # center of bin
    return times, rate


def BRER_from_brates(filename):
    data = np.loadtxt(filename)
    BR = data[:,1]
    ER = data[:,2]
    return data[:,0]-(data[1,0]-data[0,0])/2, ER, BR

def rate_from_ratefile(filename):
    data = np.loadtxt(filename)
    return data[:,0]-(data[1,0]-data[0,0])/2, data[:,1]


def average_BRER_over_realizations(filenames, binsize=20.e-3, tau=16.e-3):
    """
    Computes the sample-average event and burst rates
    :param filenames:   List of raster files
    :param binsize:     Size of the discretization bins (in seconds)
    :param tau:         Burst detection threshold (in seconds)
    :return:            Time-dependent averages event and burst rates, over different seeds for the random noise
    """
    mean_ER = np.array([])
    mean_BR = np.array([])
    for i, filename in enumerate(filenames):
        times, ER, BR = BRER_from_raster(filename, binsize=binsize, tau=tau)
        if i==0:
            mean_ER = ER
            mean_BR = BR
        else:
            mean_ER += ER
            mean_BR += BR
    mean_ER /= len(filenames)
    mean_BR /= len(filenames)
    return times, mean_ER, mean_BR


def std_BRER_over_realizations(filenames, binsize=20.e-3, tau=16.e-3):
    """
    Computes the standard deviation of event and burst rates
    :param filenames:   List of raster files
    :param binsize:     Size of the discretization bins (in seconds)
    :param tau:         Burst detection threshold (in seconds)
    :return:            Time-dependent std of event and burst rates
    """
    ERs = []
    BRs = []
    for i, filename in enumerate(filenames):
        #times, ER, BR = BRER_from_raster(filename, binsize=binsize, tau=tau)
        times, ER, BR = BRER_from_brates(filename)
        ERs.append(ER)
        BRs.append(BR)

    ERs = np.array(ERs)
    BRs = np.array(BRs)
    BPs = BRs / ERs

    mean_ER = np.mean(ERs, axis=0)
    mean_BR = np.mean(BRs, axis=0)
    mean_BP = np.mean(BPs, axis=0)

    std_ER = np.std(ERs, axis=0)
    std_BR = np.std(BRs, axis=0)
    std_BP = np.std(BPs, axis=0)

    print("std = {}".format(np.mean(std_ER)))

    std = {'ER': std_ER, 'BR': std_BR, 'BP': std_BP}
    mean = {'ER': mean_ER, 'BR': mean_BR, 'BP': mean_BP}

    return times, mean, std


def average_rate_over_realizations(filenames, binsize=20.e-3):
    """
    Computes the sample-average  rates
    :param filenames:   List of raster files
    :param binsize:     Size of the discretization bins (in seconds)
    :return:            Time-dependent averages rates, over different seeds for the random noise
    """
    mean_rate = np.array([])
    for i, filename in enumerate(filenames):
        times, rate = rates_from_raster(filename, binsize=binsize)
        if i==0:
            mean_rate = rate
        else:
            mean_rate += rate[:len(mean_rate)]

    mean_rate /= len(filenames)
    return times, mean_rate


def std_rate_over_realizations(filenames, binsize=20.e-3):
    """
    Computes the standard deviation of rates
    :param filenames:   List of raster files
    :param binsize:     Size of the discretization bins (in seconds)
    :return:            Time-dependent std of event and burst rates
    """
    rates = []
    for i, filename in enumerate(filenames):
        times, rate = rate_from_ratefile(filename)
        rates.append(rate)
        rates.append(rate)
    rates = np.array(rates)

    mean_rate = np.mean(rates, axis=0)
    std_rate = np.std(rates, axis=0)

    return times, mean_rate, std_rate


def bpercurves(fn, currents):
    """
    Compute the burst probability (BP) and event rate (ER) as a function of current.
    The simulation producing the raw data consists in step increases of currents
    applied for 1 second each. The BP and ER are computed by averaging over the last
    0.5 sec of the reponse to each step.
    :param fn:          filename (contains brate)
    :param currents:    current intensities
    :returns:
        - loc_mean_bp - list of BPs
        - loc_mean_er - list of ERs
        - loc_mean_br - list of BRs
    """
    data = np.loadtxt(fn)
    times = data[:, 0]
    bp = 100. * data[:, 1] / data[:, 2]
    br = data[:, 1]
    er = data[:, 2]
    loc_mean_bp = []
    loc_mean_er = []
    loc_mean_br = []
    for i, current in enumerate(currents):
        indices = np.where(np.logical_and(times > i + 0.5, times < i + 1.))[0]
        loc_mean_bp.append(np.mean(bp[indices]))
        loc_mean_er.append(np.mean(er[indices]))
        loc_mean_br.append(np.mean(br[indices]))
    return loc_mean_bp, loc_mean_er, loc_mean_br


def bpercurves_from_raster(fn, currents):
    """
    Compute the burst probability (BP) and event rate (ER) as a function of current.
    The simulation producing the raw data consists in step increases of currents
    applied for 1 second each. The BP and ER are computed by averaging over the last
    0.5 sec of the reponse to each step.
    :param fn:          filename (contains brate)
    :param currents:    current intensities
    :returns:
        - loc_mean_bp - list of BPs
        - loc_mean_er - list of ERs
        - loc_mean_br - list of BRs
    """
    times, er, br = BRER_from_raster(fn, binsize=20.e-3, tau=16.e-3)
    bp = 100. * br /er
    loc_mean_bp = []
    loc_mean_er = []
    loc_mean_br = []
    for i, current in enumerate(currents):
        indices = np.where(np.logical_and(times > i + 0.5, times < i + 1.))[0]
        loc_mean_bp.append(np.mean(bp[indices]))
        loc_mean_er.append(np.mean(er[indices]))
        loc_mean_br.append(np.mean(br[indices]))
    return loc_mean_bp, loc_mean_er, loc_mean_br


def ficurves(fn, currents):
    """
    Compute the firing rate as a function of current.
    The simulation producing the raw data consists in step increases of currents
    applied for 1 second each. The firing rate is computed by averaging over the last
    0.5 sec of the reponse to each step.
    :param fn:          filename (contains brate)
    :param currents:    current intensities
    :returns:
        - mean_fr - list of firing rates
    """
    data = np.loadtxt(fn)
    times = data[:, 0]
    fr = data[:, 1]
    mean_fr = []
    for i, current in enumerate(currents):
        indices = np.where(np.logical_and(times > i + 0.5, times < i + 1.))[0]
        mean_fr.append(np.mean(fr[indices]))
    return mean_fr


def get_spiketrains(rasterfile, binsize=5e-3):
    """
    Construct spike trains for all neurons in raster file
    :param rasterfile: first column constains spike times; second col contain neuron ID
    :return: spike trains (dim =  N x (number of time bins) )
    """
    raster = np.loadtxt(rasterfile)
    N_bins = int(np.max(raster[:, 0]) / binsize)
    N_neurons = int(max(raster[:, 1]))+1
    spiketrains = np.zeros((N_neurons, N_bins + 1))

    for i, spiketime in enumerate(raster[:, 0]):
        spiketrains[raster[i, 1], int(spiketime/binsize)] += 1.

    return spiketrains


def pop_corr(rasterfile, binsize=5e-3, N_selected_neurons=100):
    """
    Compute the auto- and cross-covariance of the population
    """
    spiketrains = get_spiketrains(rasterfile, binsize=binsize)

    a = 0.
    c = 0.




#################################################################
#                   PLOTTING-RELATED FUNCTIONS                  #
#################################################################
def display_BRER(filenames, outfile, binsize=20.e-3, tau=16.e-3):
    """
    Display average rates and burst probability.
    :param filenames:   List of raster files
    :param binsize:     Size of the discretization bins (in seconds)
    :param tau:         Burst detection threshold (in seconds)
    """
    t, mean, std = std_BRER_over_realizations(filenames, binsize=binsize, tau=tau)

    ER = mean['ER']
    BR = mean['BR']
    BP = 100 * mean['BP']

    std_ER = std['ER']
    std_BR = std['BR']
    std_BP = 100 * std['BP']

    plt.figure(figsize=(4, 4))
    plt.subplot(311)
    plt.plot(t, BP, color=custom_colors['red'], label='BP')
    plt.fill_between(t, BP - 2 * std_BP, BP + 2 * std_BP, color=custom_colors['red'], alpha=0.5)
    plt.ylim([0., 100.])
    plt.ylabel(r'BP [\%]')

    plt.subplot(312)
    plt.plot(t, BR, color=custom_colors['orange'], label='BR')
    plt.fill_between(t, BR - 2 * std_BR, BR + 2 * std_BR, color=custom_colors['orange'], alpha=0.5)
    plt.ylim([0., 5.])
    plt.ylabel('BR [Hz]')

    plt.subplot(313)
    plt.plot(t, ER, color=custom_colors['blue'], label='ER')
    plt.fill_between(t, ER - 2 * std_ER, ER + 2 * std_ER, color=custom_colors['blue'], alpha=0.5)
    plt.ylim([0., 13.])
    plt.xlabel('Time [s]')
    plt.ylabel('ER [Hz]')
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


def display_BRER_with_inputs(filenames, outfile, inputs, binsize=20.e-3, tau=16.e-3):
    """
     Display average rates, burst probability together with the dendritic and somatic inputs
     :param filenames:   List of raster files
     :params outfile:    Output file name
     :param inputs:      Dictionary with element 'dendrite' and 'soma' and 'times'
     :param binsize:     Size of the discretization bins (in seconds)current traces in pA
     :param tau:         Burst detection threshold (in seconds)
     """
    t, mean, std = std_BRER_over_realizations(filenames, binsize=binsize, tau=tau)

    ER = mean['ER']
    BR = mean['BR']
    BP = 100*mean['BP']

    std_ER = std['ER']
    std_BR = std['BR']
    std_BP = 100*std['BP']

    plt.figure(figsize=(4, 4))
    plt.subplot(311)
    plt.plot(t, BP, color=custom_colors['red'], label='BP')
    rescaled_input = np.min(BP) + (inputs['dendrite'] - np.min(inputs['dendrite'])) *\
                                   (np.max(BP)-np.min(BP))/(np.max(inputs['dendrite'])-np.min(inputs['dendrite']))
    plt.plot(inputs['times'], rescaled_input, '--', color=custom_colors['red'])
    plt.fill_between(t, BP - 2 * std_BP, BP + 2 * std_BP, color=custom_colors['red'], alpha=0.5)
    plt.xlim([inputs['times'][0], inputs['times'][-1]])
    plt.ylim([0., 100.])
    plt.ylabel(r'BP, $I_\mathrm{d}$ (scaled)')

    plt.subplot(312)
    plt.plot(t, BR, color=custom_colors['orange'], label='BR')
    plt.fill_between(t, BR - 2 * std_BR, BR + 2 * std_BR, color=custom_colors['orange'], alpha=0.5)
    plt.xlim([inputs['times'][0], inputs['times'][-1]])
    plt.ylim([0., 5.])
    plt.ylabel('BR [Hz]')

    plt.subplot(313)
    plt.plot(t, ER, color=custom_colors['blue'], label='ER')
    rescaled_input = np.min(ER) + (inputs['soma'] - np.min(inputs['soma'])) *\
                                   (np.max(ER)-np.min(ER))/(np.max(inputs['soma'])-np.min(inputs['soma']))
    plt.plot(inputs['times'], rescaled_input, '--', color=custom_colors['blue'])
    plt.fill_between(t, ER - 2 * std_ER, ER + 2 * std_ER, color=custom_colors['blue'], alpha=0.5)
    plt.xlim([inputs['times'][0], inputs['times'][-1]])
    plt.ylim([0., 13.])
    plt.xlabel('Time [s]')
    plt.ylabel('ER [Hz], $I_\mathrm{s}$ (scaled)')
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


def display_rates_with_inputs(filenames, outfile, inputs, encoding='event', binsize=20.e-3):
    """
     Display average rates together with the dendritic and somatic inputs
     :param filenames:   List of raster files
     :params outfile:    Output file name
     :param inputs:      Dictionary with element 'dendrite' and 'soma' and 'times'
     :param binsize:     Size of the discretization bins (in seconds)current traces in pA
     :param tau:         Burst detection threshold (in seconds)
     """

    t, rate, std_rate = std_rate_over_realizations(filenames, binsize=binsize)
    if encoding == 'event':
        current = inputs['soma']
    elif encoding == 'burst':
        current = inputs['soma']*inputs['dendrite']
    elif encoding == 'bp':
        current = inputs['dendrite']

    plt.figure(figsize=(2.5, 2.5/1.6))
    plt.plot(t, rate, color='black')
    rescaled_input = np.min(rate) + (current - np.min(current)) *\
                                   (np.max(rate)-np.min(rate))/(np.max(current)-np.min(current))
    plt.plot(inputs['times'], rescaled_input, '--', color=custom_colors['blue'])
    plt.fill_between(t, rate - 2 * std_rate, rate + 2 * std_rate, color='black', alpha=0.5)
    plt.xlim([inputs['times'][0], inputs['times'][-1]])
    plt.ylim([0., 20])
    plt.xlabel('Time [s]')
    plt.ylabel('Rate [Hz]')
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


def display_responsecurves(currents, outfile, *fns):
    """
    Plots FI curves as a function of the applied current.
    :param currents:    Currents
    :param outfile:     Figure file
    :param fns:         Files containing brate
    :return:
    """
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.5, 2.5))
    alphas = [0.3, 0.65, 1]

    for i, fn in enumerate(fns):
        bp, er, br = bpercurves(fn, currents)
        if i == 0:
            ax1.plot(currents, bp, color=custom_colors['red'], alpha=alphas[i], label='BP')
            ax1.plot(currents, er, color=custom_colors['blue'], alpha=alphas[i], label='ER')
            ax2.plot(currents, br, color=custom_colors['orange'], alpha=alphas[i], label='BR')
            ax2.plot(currents, er, color=custom_colors['blue'], alpha=alphas[i], label='ER')
        else:
            ax1.plot(currents, bp, color=custom_colors['red'], alpha=alphas[i])
            ax1.plot(currents, er, color=custom_colors['blue'], alpha=alphas[i])
            ax2.plot(currents, br, color=custom_colors['orange'], alpha=alphas[i])
            ax2.plot(currents, er, color=custom_colors['blue'], alpha=alphas[i])

    ax1.set_xlabel('$I_\mathrm{dend}$ [pA]')
    ax1.set_ylabel('BP [\%], ER [Hz]')
    ax1.set_ylim([0, 50])
    ax1.legend(loc='best')

    ax2.set_xlabel('$I_\mathrm{dend}$ [pA]')
    ax2.set_ylabel('BR [\%], ER [Hz]')
    ax2.set_ylim([0, 10])
    ax2.legend(loc='best')

    sns.despine()
    plt.tight_layout()
    #plt.savefig('BPversusIdendFromRaster_NoFB.png')
    plt.savefig(outfile)
    plt.close()

def display_ficurves(currents, outfile, *fns):
    """
    Plots FI curves as a function of the applied current.
    :param currents:    Dendritici currents
    :param outfile:     Figure file
    :param fns:         Files containing population rates
    :return:
    """
    plt.figure(figsize=(2.5*1.6, 2.5))
    alphas = [0.3, 0.65, 1]

    for i, fn in enumerate(fns):
        fr = ficurves(fn, currents)
        plt.plot(currents, fr, color='black', alpha=alphas[i])

    plt.xlabel('$I_\mathrm{dend}$ [pA]')
    plt.ylabel('Rate [Hz]')

    sns.despine()
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()



#################################################################
#                   INPUT-RELATED FUNCTIONS                     #
#################################################################
def add_step(bg, begin, length, amplitude):
    '''
    Function that adds a step of amplitude *amplitude* and duration *length*
    atop a background signal *bg*. The step starts at *begin*.
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


