from pairingprotocols import PeriodicPairing
from preeventsynapse import AdaptivePreEventSynapse
import numpy as np
import matplotlib.pyplot as plt

"""
The plasticity protocol consists in 15 sequences of 5 pairings
separated by 10 s.
"""


def simulate(freqs, eta, alpha, burst_threshold, tau_pre, starting_estimate=(5., 0.15*5.)):
    np.random.seed(1)

    # select parameters for pairing protocol
    Delta = 0.0000  # t_post - t_pre
    number_of_pairings = 5
    number_of_epochs = 15
    interepoch_interval = 10.

    # integration time step
    dt = 0.0001

    weights = np.zeros(len(freqs))

    for i, f in enumerate(freqs):

        last_estimated_ER = starting_estimate[0]
        last_estimated_BR = starting_estimate[1]

        for epoch in range(number_of_epochs):
            syn = AdaptivePreEventSynapse(eta, tau_trace=tau_pre, tau_ma=alpha, burst_def=burst_threshold,
                                          starting_estimate=(last_estimated_ER, last_estimated_BR))

            protocol = PeriodicPairing(f, Delta, pairing_number=number_of_pairings)

            T = number_of_pairings / f + 1.

            syn.pre.compute_trains(protocol.spiketimes_pre, dt, T)
            syn.post_ma.compute_trains(protocol.spiketimes_post, dt, T)

            for t in range(int(T/dt)):
                syn.evolve(t, dt)
            weights[i] += syn.weight/dt
            #  division by dt is because event and burst trains are not divided
            #  by dt in the rule (see preeventsynapse.py)

            last_estimated_ER = syn.post_ma.moving_average['event']
            last_estimated_ER = last_estimated_ER*np.exp(-interepoch_interval/alpha)

            last_estimated_BR = syn.post_ma.moving_average['burst']
            last_estimated_BR = last_estimated_BR*np.exp(-interepoch_interval/alpha)

    return freqs, weights


if __name__ == '__main__':
    import sys
    import seaborn as sns
    sys.path.append('../')
    plt.style.use('../thesis_mplrc.dms')

    freqs = np.arange(10, 100., 10)  # frequency in Hz
    alpha = 15.
    tau_pre = 0.1
    f, w = simulate(freqs, 0.1, alpha, 0.016, tau_pre)

    plt.figure(figsize=(3, 3 / 1.6))
    plt.plot(f, w, color='black', label='simul.')
    plt.plot(f, np.zeros(f.shape), '--k', lw=0.5)
    sns.despine()
    plt.xlabel('Pairing frequency [Hz]')
    plt.ylabel('$\Delta W$')
    plt.tight_layout()
    plt.savefig('../../results/learning-rule/Periodic/DeltaW_AdaptiveLearningRule_Alpha'+str(alpha)+'_Taupre'+str(tau_pre)+'.pdf')
    plt.close()