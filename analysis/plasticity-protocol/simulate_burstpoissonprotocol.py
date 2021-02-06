from pairingprotocols import PoissonBurstProtocol
from preeventsynapse import AdaptivePreEventSynapse
import numpy as np
import matplotlib.pyplot as plt

"""
Description:
The plasticity protocol consists in pre and postsynaptic Poisson event-burst
processes. Namely, a neuron fires an event, and a burst is then evoked with probability p.
"""


def simulate(bps, eta, duration, alpha, burst_threshold, tau_pre, event_rate):
    """
    Compute the weight change for the burst-poisson protocol, which consists in pre and postsynaptic Poisson event-burst
processes. Namely, a neuron fires an event, and a burst is then evoked with probability p.
    :param bps:             burst probabilities of the protocol
    :param eta:             learning rate
    :param duration:        duration of the protocol (s)
    :param alpha:           time constant of exponential moving average
    :param burst_threshold: threshold for burst detection
    :param tau_pre:         decay time constant of presynaptic trace
    :param event_rate:      rate of event generation
    :return:                weights
    """
    # Set random seed for reproducibility
    np.random.seed(2)

    # Parameters
    nb_reals = 20       # number of realizations of the random pairing
    ER = event_rate     # event rate
    bp_protocol = 0.2   # initial BP estimate
    dt = 0.001          # integration time step

    # Main simulation
    weights = np.zeros((nb_reals, len(bps)))
    for n in range(nb_reals):
        if (n+1) % 10 == 0:
            print('realization #{}'.format(n+1))
        for i, p in enumerate(bps):
            syn = AdaptivePreEventSynapse(eta,
                                          tau_trace=tau_pre,
                                          tau_ma=alpha,
                                          burst_def=burst_threshold,
                                          starting_estimate=(ER, bp_protocol*ER))

            protocol = PoissonBurstProtocol(p=p,
                                            duration=duration,
                                            rate_event=ER/(1 - ER*0.020))
            # Note: rate_event is equal to ER/(1 - ER*abs_ref_event) to yield the event rate = ER in
            # the protocol. This is not essential a priori.

            syn.pre.compute_trains(protocol.spiketimes_pre, dt, duration)
            syn.post_ma.compute_trains(protocol.spiketimes_post, dt, duration)

            for t in range(int(duration/dt)):
                syn.evolve(t, dt)
            weights[n, i] += syn.weight/dt
            #  division by dt is because event and burst trains are not divided
            #  by dt in the rule (see preeventsynapse.py)

    return np.mean(weights, axis=0), np.std(weights, axis=0)


if __name__ == '__main__':
    import sys
    sys.path.append('../')
    plt.style.use('../thesis_mplrc.dms')
    import seaborn as sns

    # parameters
    tau = np.arange(0., 0.050, 0.0001)
    rates = np.linspace(1, 50, 10)
    duration = 60.
    bps = np.arange(0., 0.5, 0.05)
    alpha = 30.

    # simulate
    w = simulate(bps, 0.1, duration, alpha, 0.016, 0.020)

    # plotting
    plt.figure(figsize=(3, 3 / 1.6))
    plt.plot(100*bps, w, 'k')
    plt.plot(100*bps, np.zeros(bps.shape), 'k--', lw=1)
    sns.despine()
    plt.xlabel('Initial BP [\%]')
    plt.ylabel('$\Delta W$')
    plt.tight_layout()
    plt.savefig('../../results/learning-rule/Poisson/DeltaW_BurstPoisson_Alpha' + str(alpha) + '.pdf')
    plt.close()