from pairingprotocols import PeriodicPairing
from preeventsynapse import AdaptivePreEventSynapse
import numpy as np
import matplotlib.pyplot as plt

"""
Description:
The plasticity protocol consists in 15 bursts of 5 pairings, assuming a long
time in between bursts. This way, we can simply consider 5 pairings
and multiply the weight change by 15 at the end (since there is no weight decay).
"""

def simulate(freqs, alpha):
    np.random.seed(1)

    # select parameters for pairing protocol
    Delta = 0.0000  # t_post - t_pre
    number_of_pairings = 5
    number_of_epochs = 15
    interepoch_interval = 10.

    # create synapse
    starting_estimate = (5, 0.15*5.)

    # integration time step
    dt = 0.0001

    weights = np.zeros(len(freqs))

    for i, f in enumerate(freqs):

        last_estimated_ER = starting_estimate[0]
        last_estimated_BR = starting_estimate[1]

        for epoch in range(number_of_epochs):
            syn = AdaptivePreEventSynapse(0.1, tau_trace=0.020, tau_ma=alpha, burst_def=0.016,
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
            #print("estimated BP = {}".format(last_estimated_BR/last_estimated_ER))

    return freqs, weights


def analytical(eta, A_plus, A_minus, burst_threshold):
    freqs = np.arange(2, 150., 2.)  # frequency in Hz
    Delta = 0.000  # t_post - t_pre
    tau_pre = 0.005
    assert(np.all(np.floor(freqs*Delta)) == 0)

    a = np.exp(-1./(freqs*tau_pre))

    weights = (freqs < 1/burst_threshold)*15.*eta*\
              (A_minus/A_plus)*( np.exp(-Delta/tau_pre)/(a-1.) )*( a*(a**5-1)/(a-1.) - 5. ) + (freqs > 1 / burst_threshold)*15*eta*(a + A_minus/A_plus)

    return freqs, weights


if __name__ == '__main__':
    import sys
    sys.path.append('../')
    plt.style.use('../thesis_mplrc.dms')
    import seaborn as sns

    freqs = np.arange(2, 80., 2.)  # frequency in Hz
    alpha = 10.
    f, w = simulate(freqs, alpha)
    #f_anal, w_anal = analytical(0.1, 1., -0.15, 0.016)

    plt.figure(figsize=(3, 3 / 1.6))
    plt.plot(f,w, color='black', label='simul.')
    #plt.plot(f_anal, w_anal, label='anal.')
    plt.plot(f, np.zeros(f.shape), '--k', lw=0.5)
    sns.despine()
    plt.xlabel('Pairing frequency [Hz]')
    plt.ylabel('$\Delta W$')
    #plt.legend()
    plt.tight_layout()
    plt.savefig('../../results/learning-rule/Periodic/DeltaW_AdaptiveLearningRule_Alpha' + str(alpha) + '.pdf')
    plt.close()