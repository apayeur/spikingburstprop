from pairingprotocols import PoissonBurstProtocol
from preeventsynapse import AdaptivePreEventSynapse
import numpy as np
import matplotlib.pyplot as plt

"""
Description:
The plasticity protocol consists in pre- and post-synaptic Poisson event-burst
processes. Namely, a neuron fires an event, and a burst is then evoked with probability p.
Here, the plasticity rule is adaptive, i.e., the postsynaptic neuron
estimates its burst probability.
"""


def simulate(bps, eta, duration, alpha, burst_threshold, tau_pre):
    """
    Compute the weight change for the burst-poisson protocol for
    different initial estimated
    :param bps:             Pre-pairing burst probabilities (in [0,1])
    :param eta:             Learning rate
    :param duration:        Duration of the protocol (s)
    :param alpha:           Time constant of exponential moving average
    :param burst_threshold: Threshold for burst detection
    :param tau_pre:         Decay time constant of presynaptic trace
    :return:                Weights
    """
    # set random seed for reproducibility
    np.random.seed(7)

    # parameters
    nb_reals = 20       # number of realizations of the random pairing
    ER = 5.             # event rate
    bp_protocol = 0.2  # Burst probability used to generate bursts in protocol
    dt = 0.001          # integration time step

    weights = np.zeros(bps.shape)

    for n in range(nb_reals):
        if (n+1) % 2 == 0:
            print('realization #{}'.format(n+1))
        for i, p in enumerate(bps):
            syn = AdaptivePreEventSynapse(eta, tau_trace=tau_pre, tau_ma=alpha, burst_def=burst_threshold,
                                          starting_estimate=(ER, p*ER))

            protocol = PoissonBurstProtocol(p=bp_protocol, duration=duration, rate_event=ER/(1 - ER*0.020))
            # Note: rate_event is equal to ER/(1 - ER*abs_ref_event) to yield the event rate = ER in
            # the protocol. This is not essential a priori.

            syn.pre.compute_trains(protocol.spiketimes_pre, dt, duration)
            syn.post_ma.compute_trains(protocol.spiketimes_post, dt, duration)

            for t in range(int(duration/dt)):
                syn.evolve(t, dt)
            weights[i] += syn.weight/dt
            #  division by dt is because event and burst trains are not divided
            #  by dt in the rule (see preeventsynapse.py)
    weights /= nb_reals

    return weights


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
    #bp_anal, w_anal = analytical(0.1, 400, 5., 1, -0.15, 0.016, 0.020, 0.020)
    #print(burst_fraction(np.linspace(1., 30., 10), 0.016))

    # plotting
    plt.figure(figsize=(3, 3 / 1.6))
    plt.plot(100*bps, w, 'k', label='simul.')
    #plt.plot(bp_anal, w_anal, label='anal.', color='grey')
    plt.plot(100*bps, np.zeros(bps.shape), 'k--', lw=1)
    #plt.legend()
    sns.despine()
    plt.xlabel('Initial BP [\%]')
    plt.ylabel('$\Delta W$')
    plt.tight_layout()
    plt.savefig('../../results/learning-rule/Poisson/DeltaW_AdaptiveLearningRule_BurstPoisson_Alpha' + str(alpha) + '.pdf')
    plt.close()