from pairingprotocols import PeriodicPairing
from preeventsynapse import PreEventSynapse
import numpy as np
import matplotlib.pyplot as plt

"""
Description:
The plasticity protocol consists in 15 bursts of 5 pairings, assuming a long
time in between bursts. This way, we can simply consider 5 pairings
and multiply the weight change by 15 at the end (since there is no weight decay).
"""

def simulate():
    np.random.seed(1)

    # select parameters for pairing protocol
    freqs = np.arange(2, 80., 2.)  # frequency in Hz
    Delta = 0.0000  # t_post - t_pre
    number_of_pairings = 5

    # create synapse
    syn = PreEventSynapse(0.1, tau_trace=0.009)

    # integration time step
    dt = 0.0001

    weights = []
    for f in freqs:
        syn.reset()
        protocol = PeriodicPairing(f, Delta, pairing_number=number_of_pairings)
        #protocol.display('../../results/learning-rule/protocol_f'+str(f)[:2]
        #                 +'_Delta'+ str(Delta*1000)[:2]+'.pdf')
        T = number_of_pairings / f + 1.
        syn.pre.compute_trains(protocol.spiketimes_pre, dt, T)
        syn.post.compute_trains(protocol.spiketimes_post, dt, T)

        for t in range(int(T/dt)):
            syn.evolve(t, dt)
        weights.append(15*syn.weight/dt)
        #  division by dt is because event and burst trains are not divided
        #  by dt in the rule (see preeventsynapse.py)
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
    #f, w = simulate()
    f_anal, w_anal = analytical(0.1, 1., -0.15, 0.016)
    #plt.plot(f,w, label='simul.')
    plt.plot(f_anal, w_anal, label='anal.')
    plt.plot(f_anal, np.zeros(f_anal.shape), '--k', lw=0.5)

    plt.legend()
    plt.tight_layout()
    plt.show()