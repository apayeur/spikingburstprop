from pairingprotocols import PoissonBurstProtocol
from preeventsynapse import PreEventSynapse
import numpy as np
import matplotlib.pyplot as plt

"""
Description:
The plasticity protocol consists in pre- and post-synaptic Poisson event-burst
processes. Namely, a neuron fires an event, and a burst is then evoked with probability p.
"""


def simulate():
    np.random.seed(7)

    # select parameters for pairing protocol
    burst_probs = np.arange(0., 0.5, 0.05)
    duration = 400.  # in seconds

    # create synapse
    learning_rate = 0.1
    syn = PreEventSynapse(learning_rate, tau_trace=20.e-3)

    # integration time step
    dt = 0.001

    weights = []
    for p in burst_probs:
        syn.reset()
        protocol = PoissonBurstProtocol(p=p, duration=duration, rate_event=5.)

        # uncomment the following line to output the spiking structure of the protocol
        #protocol.display('../../results/learning-rule/Poisson/protocol_p'+str(p)[:4]+'.pdf')
        syn.pre.compute_trains(protocol.spiketimes_pre, dt, duration)
        syn.post.compute_trains(protocol.spiketimes_post, dt, duration)

        for t in range(int(duration/dt)):
            syn.evolve(t, dt)
        print(syn)
        weights.append(syn.weight/dt)
        #  division by dt is because event and burst trains are not divided
        #  by dt in the rule (see preeventsynapse.py)

    return burst_probs, weights

def analytical(eta, duration, rate, A_plus, A_minus, burst_threshold, tau_pre, tau_ref_E):
    burst_probs = np.arange(0., 0.5, 0.05)
    return burst_probs, duration*eta*tau_pre*(burst_probs + A_minus/A_plus)*(rate/(1.+rate*tau_ref_E))**2


if __name__ == '__main__':
    bp, w = simulate()
    bp_anal, w_anal = analytical(0.1, 400, 5., 1, -0.15, 0.016, 0.020, 0.020)
    #print(burst_fraction(np.linspace(1., 30., 10), 0.016))
    plt.plot(bp,w, label='simul.')
    plt.plot(bp_anal, w_anal, label='anal.', color='grey')
    plt.plot(bp, np.zeros(bp.shape), 'k--', lw=1)
    plt.legend()
    plt.tight_layout()
    plt.show()