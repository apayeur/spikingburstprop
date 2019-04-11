import matplotlib.pyplot as plt
from custom_colors import custom_colors
import seaborn as sns
import numpy as np

plt.style.use('/Users/alexandrepayeur/Documents/Recherches/Burstprop/Codes/thesis_mplrc.dms')


class PairingProtocol(object):
    """
    A generic pairing protocol for a plasticity experiment.

    Attributes:
        spiketimes_pre: list of forced presynaptic spike times (float, in seconds)
        spiketimes_post: list of forced postsynaptic spike times (float, in seconds)
    """

    def __init__(self):
        """Return a pairing protocol with *pairing_number* pre-post spike pairs
        and empty lists of presynaptic spike times *spiketimes_pre* and postsynaptic spike
        times *spiketimes_post*.
        """
        self.spiketimes_pre = []
        self.spiketimes_post = []

    def display(self, filename='../../results/fig3-learning-rule/protocol.pdf'):
        """
        Plot pairing protocol
        """
        try:
            assert(len(self.spiketimes_pre)>0 and len(self.spiketimes_post)>0) # check if spiketimes list are non-empty
            plt.figure(figsize=(5,2))
            plt.eventplot(self.spiketimes_pre, colors=[custom_colors['red']], lineoffsets=1, linelengths=2, lw=0.5, label='pre')
            plt.eventplot(self.spiketimes_post, colors=[custom_colors['blue']], lineoffsets=2, linelengths=2, lw=0.5, label='post')
            plt.plot([0.5, 0.510], [3.5, 3.5], lw=0.5)
            plt.xlim([0,min(5, max(self.spiketimes_post[-1], self.spiketimes_pre[-1]))])
            #plt.xticks([])
            #plt.yticks([])
            sns.despine(left=True, bottom=True)
            plt.savefig(filename)
            plt.close()
        except:
            print('Empty spiketimes list. Nothing to plot.')


class PeriodicPairing(PairingProtocol):
    """
       Periodic pairing protocol for a plasticity experiment.

       Attributes:
           pairing_number: number of pairings (integer)
           frequency: pairing frequency (float, in Hz)
           Deltat: time difference between pre and post spiks (float, in s)
       """

    def __init__(self, frequency, Deltat, pairing_number=60):
        """Return a pairing protocol with periodic pre-post pairing.
        These *pairing_number* pairings occur with frequency *frequency*
        and with a time difference *Deltat*."""
        super().__init__()
        self.pairing_number = pairing_number
        self.frequency = frequency
        self.Deltat = Deltat
        self.get_spiketimes()

    def get_spiketimes(self):
        """
        Compute the pre and post spike times and put them in lists.
        """
        delay = 0.010  # delay before first pairing

        if len(self.spiketimes_pre)>0 and len(self.spiketimes_post)>0:
            self.spiketimes_pre = []
            self.spiketimes_post = []

        if self.Deltat > 0:
            for i in range(self.pairing_number):
                self.spiketimes_pre.append(delay + i/self.frequency)
                self.spiketimes_post.append(delay + i/self.frequency + self.Deltat)
        else:
            for i in range(self.pairing_number):
                self.spiketimes_post.append(delay + i/self.frequency)
                self.spiketimes_pre.append(delay + i/self.frequency - self.Deltat)


class PoissonBurstProtocol(PairingProtocol):
    """
       Pre- and post-synaptic neurons fire events with Poisson rate rate_event.
       If an event occurs, then a burst is evoked with probability p.

       Attributes:
           duration : time duration of protocol (float, in seconds)
           rate_event: rate of pre- and post-synaptic events (float, in Hz)
           p: burst propbability (float)
    """

    def __init__(self, duration=10., p=0.15, rate_event=5.):
        """Return a pairing protocol with Poisson events at rate *rate_event*
        converted into bursts with probability *p*
        """
        super().__init__()
        self.duration = duration
        self.p = p
        self.rate_event = rate_event
        self.abs_ref_per = 0.004  # absolute refractory due to spike shape
        self.abs_ref_per_event = 0.020  # absolute refractory for events
        self.get_spiketimes()

    def get_spiketimes(self):
        """
        Compute the pre and post spike times and put them in lists.
        """
        if len(self.spiketimes_pre) > 0 and len(self.spiketimes_post) > 0:
            self.spiketimes_pre = []
            self.spiketimes_post = []

        last_spike = 0.
        while True:
            interval = self.abs_ref_per_event + np.random.exponential(1./self.rate_event)
            last_spike = last_spike + interval
            if last_spike > self.duration:
                break
            self.spiketimes_pre.append(last_spike)

            r = np.random.uniform(0., 1.)
            if r < self.p:  # check if burst
                interval = self.abs_ref_per + np.random.uniform(0., 0.010)
                last_spike = last_spike + interval
                if last_spike > self.duration:
                    break
                self.spiketimes_pre.append(last_spike)

        last_spike = 0.
        while True:
            interval = self.abs_ref_per_event + np.random.exponential(1. / self.rate_event)
            last_spike = last_spike + interval
            if last_spike > self.duration:
                break
            self.spiketimes_post.append(last_spike)

            r = np.random.uniform(0., 1.)
            if r < self.p:  # check if burst
                interval = self.abs_ref_per + np.random.uniform(0., 0.010)
                last_spike = last_spike + interval
                if last_spike > self.duration:
                    break
                self.spiketimes_post.append(last_spike)


class RandomProtocol(PairingProtocol):
    """
       Pre- and post-synaptic neurons firing with Poisson rate.
    """

    def __init__(self, spike_number=75, p=0.15, rate=5.):
        """Return a pairing protocol with Poisson events at rate *rate*.
        """
        super().__init__()
        self.spike_number = spike_number
        self.p = p
        self.rate = rate
        self.get_spiketimes()

    def get_spiketimes(self):
        """
        Compute the pre and post spike times and put them in lists.
        """
        if len(self.spiketimes_pre) > 0 and len(self.spiketimes_post) > 0:
            self.spiketimes_pre = []
            self.spiketimes_post = []

        isis = np.random.exponential(1./self.rate, self.spike_number)
        self.spiketimes_pre = list(isis.cumsum())
        isis = np.random.exponential(1./self.rate, self.spike_number)
        self.spiketimes_post = list(isis.cumsum())


class SjostromProtocol(PairingProtocol):
    """
       15 bursts of 5 pre- and post-synaptic neurons spikes with Poisson rate.
    """

    def __init__(self, duration=10., p=0.15, rate=5.):
        """Return a pairing protocol with Poisson events at rate *rate*.
        """
        super().__init__()
        self.duration = duration
        self.p = p
        self.rate = rate
        self.get_spiketimes()

    def get_spiketimes(self):
        """
        Compute the pre and post spike times and put them in lists.
        """
        if len(self.spiketimes_pre) > 0 and len(self.spiketimes_post) > 0:
            self.spiketimes_pre = []
            self.spiketimes_post = []

        last_spike = 0.
        while True:
            interval = np.random.exponential(1./self.rate)
            last_spike = last_spike + interval
            if last_spike > self.duration:
                break
            self.spiketimes_pre.append(last_spike)

        last_spike = 0.
        while True:
            interval = np.random.exponential(1. / self.rate)
            last_spike = last_spike + interval
            if last_spike > self.duration:
                break
            self.spiketimes_post.append(last_spike)
