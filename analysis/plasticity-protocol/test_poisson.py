import numpy as np
import matplotlib.pyplot as plt

# parameters
synaptic_params = dict()
synaptic_params['presynaptic trace time const'] = 0.050
synaptic_params['moving average time const']    = 15.
synaptic_params['learning rate']                = 0.1
synaptic_params['burst detection threshold']    = 0.016

protocol_params = dict()
protocol_params['time step']              = 0.001
protocol_params['duration']               = 60.
protocol_params['number of realizations'] = 10
protocol_params['rates']                  = np.linspace(1., 50., 10)
protocol_params['Pbar0']                  = 0.5


def poisson_spike_times(rate, duration):
    N = int(4*rate*duration)
    while True:
        presyn_times = np.cumsum(np.random.exponential(1./rate, N))
        postsyn_times = np.cumsum(np.random.exponential(1./rate, N))
        if min(presyn_times[-1], postsyn_times[-1]) > duration:
            break
        else:
            N += rate*duration
    return presyn_times[np.where(presyn_times<duration)[0]], postsyn_times[np.where(postsyn_times<duration)[0]]


def extract_event_times(spiketimes):
    burst_state = -1
    eventtimes = []
    bursttimes = []
    eventtimes.append(spiketimes[0])
    for i in range(1, len(spiketimes)):
        if spiketimes[i] - spiketimes[i - 1] > synaptic_params['burst detection threshold']:
            burst_state = -1.
            eventtimes.append(spiketimes[i])
        else:
            if burst_state < 0:
                burst_state = 1
                bursttimes.append(spiketimes[i])
    return np.array(eventtimes), np.array(bursttimes)


def compute_pbar(eventtimes, bursttimes):
    et, bt = construct_event_trains(eventtimes, bursttimes)
    mean_e = len(eventtimes)/protocol_params['duration']
    t = np.arange(0, 10*synaptic_params['moving average time const'], protocol_params['time step'])
    alpha_filter = (1./synaptic_params['moving average time const'])*np.exp(-t/synaptic_params['moving average time const'])
    Ebar = np.convolve(et, alpha_filter, 'full')[:int(protocol_params['duration']/protocol_params['time step'])]
    Bbar = np.convolve(bt, alpha_filter, 'full')[:int(protocol_params['duration']/protocol_params['time step'])]
    t = np.arange(0, protocol_params['duration'], protocol_params['time step'])
    Pbar = (protocol_params['Pbar0']*mean_e*np.exp(-t/synaptic_params['moving average time const']) + Bbar)/(mean_e*np.exp(-t/synaptic_params['moving average time const']) + Ebar)
    #plt.plot(t, Ebar[:int(protocol_params['duration']/protocol_params['time step'])])
    #plt.plot(t, Bbar[:int(protocol_params['duration']/protocol_params['time step'])])

    #plt.plot(t, Bbar/Ebar)
    #plt.show()
    return Pbar


def construct_event_trains(eventtimes, bursttimes):
    event_train = np.zeros(int(protocol_params['duration'] / protocol_params['time step']))
    burst_train = np.zeros(int(protocol_params['duration'] / protocol_params['time step']))
    event_train[(eventtimes/protocol_params['time step']).astype(int)] = 1.
    burst_train[(bursttimes/protocol_params['time step']).astype(int)] = 1.
    return event_train, burst_train


def construct_event_train(eventtimes):
    event_train = np.zeros(int(protocol_params['duration'] / protocol_params['time step']))
    event_train[(eventtimes/protocol_params['time step']).astype(int)] = 1.
    return event_train


def compute_Etilde(eventtimes):
    et = construct_event_train(eventtimes)
    t = np.arange(0, 10.*synaptic_params['presynaptic trace time const'], protocol_params['time step'])
    tilde_filter = np.exp(-t / synaptic_params['presynaptic trace time const'])
    return np.convolve(et, tilde_filter, 'full')[:int(protocol_params['duration'] / protocol_params['time step'])]


# simulation
weight_changes = np.zeros(protocol_params['rates'].shape)
BR = np.zeros(protocol_params['rates'].shape)
ER = np.zeros(protocol_params['rates'].shape)
SR = np.zeros(protocol_params['rates'].shape)
for n in range(protocol_params['number of realizations']):
    print("computing realization {}....".format(n+1))
    for i, r in enumerate(protocol_params['rates']):
        preST, postST = poisson_spike_times(r, protocol_params['duration'])
        #print(np.diff(postST)[-10:])
        #print(postST[-10:])

        #print('mean rate = {}'.format(len(postST)/protocol_params['duration']))
        preET, preBT = extract_event_times(preST)

        postET, postBT = extract_event_times(postST)
        #ER[i] = len(postET)/protocol_params['duration']
        #BR[i] = len(postBT)/protocol_params['duration']
        #SR[i] = len(postST)/protocol_params['duration']
        Pbar = compute_pbar(postET, postBT)
        Etilde = compute_Etilde(preET)
        E, B = construct_event_trains(postET, postBT)

        dwdt = synaptic_params['learning rate']*(B - Pbar*E)*Etilde
        Deltaw = np.sum(dwdt)
        weight_changes[i] += Deltaw
weight_changes /= protocol_params['number of realizations']


if __name__ == '__main__':
    import sys
    sys.path.append('../')
    plt.style.use('../thesis_mplrc.dms')
    import seaborn as sns
    # plotting
    plt.figure(1, figsize=(4,3))
    plt.plot(protocol_params['rates'], weight_changes)
    plt.savefig('../../results/learning-rule/test_poisson.pdf')
    plt.close()
    #plt.figure(1, figsize=(4,3))
    #plt.plot(rates, ER, label='ER')
    #plt.plot(rates, BR, label='BR')
    #plt.plot(rates, SR, label='SR')
    #plt.plot(rates, BR/ER*100)
    #plt.legend()
#plt.show()