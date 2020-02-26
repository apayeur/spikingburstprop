import numpy as np
from scipy import signal


class NeuronState(object):
    def __init__(self, burst_def=0.016, tau_trace=0.020):
        self.burst_def = burst_def  # timescale to distinguish bursts from events (s)
        self.tau_trace = tau_trace  # time constant of exponential decay of traces
        self.train = {'event': None, 'burst': None, 'spike': None}
        self.trace = {'event': 0., 'burst': 0, 'spike': 0.}

    def compute_trains(self, spiketimes, timestep, T):
        assert (len(spiketimes) > 0)
        assert (T > spiketimes[-1])

        self.train['spike'] = np.zeros(int(T / timestep)+1)
        self.train['event'] = np.zeros(int(T / timestep)+1)
        self.train['burst'] = np.zeros(int(T / timestep)+1)

        self.train['spike'][int(spiketimes[0] / timestep)] += 1.
        self.train['event'][int(spiketimes[0] / timestep)] += 1.
        burst_state = -1

        for s in range(1, len(spiketimes)):
            self.train['spike'][int(spiketimes[s] / timestep)] += 1.
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

    def evolve_traces_except_event_trace(self, timeindex, timestep):
        self.trace['spike'] += -self.trace['spike']*(timestep/self.tau_trace) + self.train['spike'][timeindex]
        self.trace['burst'] += -self.trace['burst']*(timestep/self.tau_trace) + self.train['burst'][timeindex]

    def evolve_event_trace(self, timeindex, timestep):
        self.trace['event'] += -self.trace['event']*(timestep/self.tau_trace) + self.train['event'][timeindex]

    def reset(self):
        self.train = {'event': None, 'burst': None, 'spike': None}
        self.trace = {'event': 0., 'burst': 0, 'spike': 0.}

    def compute_eb_corr(self, dt):
        assert len(self.train['burst']) and len(self.train['event'])
        ebcorr_be = signal.correlate(np.array(self.train['burst'])/dt,
                                  np.array(self.train['event'])/dt, mode='full', method='fft')
        ebcorr_eb = signal.correlate(np.array(self.train['event'])/dt,
                                  np.array(self.train['burst'])/dt, mode='full', method='fft')
        L = len(self.train['burst'])
        index_for_zero_lag = L - 1
        k = np.arange(0, 2 * L - 1)
        k = k - index_for_zero_lag
        return k, ebcorr_be/(L - abs(k)), ebcorr_eb/(L - abs(k))

    def compute_eb_cov(self, dt):
        k, corr_be, corr_eb = self.compute_eb_corr(dt)
        ebcov_be = corr_be - np.mean(self.train['burst'])*np.mean(self.train['event'])/dt**2
        ebcov_eb = corr_eb - np.mean(self.train['burst'])*np.mean(self.train['event'])/dt**2
        return k, ebcov_be, ebcov_eb


class NeuronStateWithMovingAverage(NeuronState):
    def __init__(self, burst_def=0.016, tau_trace=0.020, tau_ma=20., starting_estimate=(5, 0.15*5)):
        super().__init__(burst_def=burst_def, tau_trace=tau_trace)
        self.tau_moving_average = tau_ma
        self.starting_estimate = starting_estimate
        self.moving_average = {'event': starting_estimate[0], 'burst': starting_estimate[1]}

    def evolve_moving_averages(self, timeindex, timestep):
        self.moving_average['event'] += -self.moving_average['event'] * (timestep / self.tau_moving_average) \
                                        + self.train['event'][timeindex] / self.tau_moving_average
        self.moving_average['burst'] += -self.moving_average['burst'] * (timestep / self.tau_moving_average) \
                                        + self.train['burst'][timeindex] / self.tau_moving_average

    def reset(self):
        super().reset()
        self.moving_average = {'event': self.starting_estimate[0], 'burst': self.starting_estimate[1]}

    def set_starting_estimate(self, s_est):
        self.starting_estimate = s_est
        self.moving_average = {'event': s_est[0], 'burst': s_est[1]}
