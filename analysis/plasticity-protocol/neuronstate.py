import numpy as np


class NeuronState(object):
    def __init__(self, burst_def=0.016, tau_trace=0.020):
        self.burst_def = burst_def  # timescale to distinguish bursts from events (s)
        self.tau_trace = tau_trace  # time constant of exponential decay of traces
        self.train = {'event': None, 'burst': None, 'spike': None}
        self.trace = {'event': 0., 'burst': 0, 'spike': 0.}

    def compute_trains(self, spiketimes, timestep, T):
        assert (T > spiketimes[-1])
        self.train['spike'] = np.zeros(int(T / timestep))
        self.train['event'] = np.zeros(int(T / timestep))
        self.train['burst'] = np.zeros(int(T / timestep))

        for spiketime in spiketimes:
            self.train['spike'][int(spiketime / timestep)] += 1.

        if len(spiketimes) > 0:
            self.train['event'][int(spiketimes[0] / timestep)] += 1.
            burst_state = -1
            for s in range(1, len(spiketimes)):
                if spiketimes[s] - spiketimes[s - 1] > self.burst_def:
                    burst_state = -1.
                    self.train['event'][int(spiketimes[s] / timestep)] += 1.
                else:
                    if burst_state < 0:
                        burst_state = 1
                        self.train['burst'][int(spiketimes[s] / timestep)] += 1.

    def evolve_traces(self, timeindex, timestep):
        self.trace['spike'] += -self.trace['spike']*(timestep/self.tau_trace) + self.train['spike'][timeindex]
        self.trace['event'] += -self.trace['event']*(timestep/self.tau_trace) + self.train['event'][timeindex]
        self.trace['burst'] += -self.trace['burst']*(timestep/self.tau_trace) + self.train['burst'][timeindex]

    def reset(self):
        self.train = {'event': None, 'burst': None, 'spike': None}
        self.trace = {'event': 0., 'burst': 0, 'spike': 0.}