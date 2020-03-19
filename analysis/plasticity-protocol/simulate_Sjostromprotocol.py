from pairingprotocols import SjostromProtocol
from preeventsynapse import AdaptivePreEventSynapse
import numpy as np
import matplotlib.pyplot as plt

"""
Description:
This plasticity protocol is described in Sjostrom et al. Neuron (2001).
It consists in 15 bursts of 5 pre and postsynaptic neurons spikes with Gaussian jitter.
Here, the plasticity rule is adaptive, i.e., the postsynaptic neuron
estimates its burst probability.
"""


def simulate(freqs, eta, alpha, burst_threshold, tau_pre, starting_estimate=(5., 0.15*5.)):
    np.random.seed(2)

    # select parameters for pairing protocol
    nb_reals = 20
    nb_pairings = 60

    # create synapse
    syn = AdaptivePreEventSynapse(eta,
                                  burst_def=burst_threshold,
                                  tau_trace=tau_pre,
                                  tau_ma=alpha,
                                  starting_estimate=starting_estimate)

    # integration time step
    dt = 0.001

    weights = np.zeros(freqs.shape)
    for n in range(nb_reals):
        print('simulating realization #{}....'.format(n+1))
        for idx, f in enumerate(freqs):
            syn.reset()
            protocol = SjostromProtocol(frequency=f, pairing_number=nb_pairings)
            #protocol.display(filename='../../results/SjostromProtocol.pdf')
            duration = protocol.get_delay() + (nb_pairings - 1)/f + 10*(int(nb_pairings/5) - 1) + 1
            syn.pre.compute_trains(protocol.spiketimes_pre, dt, duration)
            syn.post_ma.compute_trains(protocol.spiketimes_post, dt, duration)

            for x in range(int(duration/dt)):
                syn.evolve(x, dt)
            weights[idx] += syn.weight/dt   # division by dt is because event and burst trains are not divided
                                            # by dt in the rule (see preeventsynapse.py)
    weights /= nb_reals

    return freqs, weights


if __name__ == '__main__':
    import sys
    sys.path.append('../')
    plt.style.use('../thesis_mplrc.dms')
    import seaborn as sns

    freqs = np.linspace(1., 50, 4)
    pbar0 = 0.5
    alpha = 30
    tau_pre = 0.1

    for pbar0 in [0.25, 0.5, 0.75]:
        for tau_pre in [0.05, 0.2, 0.5]:
            for alpha in [5, 10, 30]:
                r, w = simulate(freqs, 0.1, alpha, 0.016, tau_pre, starting_estimate=(2., pbar0*2.))
                plt.figure(figsize=(3, 3/1.6))
                plt.plot(r, w, color='black')
                sns.despine()
                plt.plot(r, np.zeros(r.shape), 'k--', lw=1)
                plt.xlabel('Frequency [Hz]')
                plt.ylabel('$\Delta W$')
                plt.tight_layout()
                plt.savefig('../../results/learning-rule/DeltaW_AdaptiveLearningRule_SjostromProtocol_Alpha'
                            +str(alpha)+'_Taupre'+str(tau_pre)+'_pbar0'+str(pbar0)+'.pdf')
                plt.close()