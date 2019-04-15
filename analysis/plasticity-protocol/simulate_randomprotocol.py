from pairingprotocols import RandomProtocol
from preeventsynapse import PreEventSynapse
import numpy as np
import matplotlib.pyplot as plt

"""
Description:
The plasticity protocol consists in pre- and post-synaptic Poisson processes.
"""

def simulate():
    np.random.seed(1)

    # select parameters for pairing protocol
    rates = np.linspace(1., 20., 10)
    duration = 400 # seconds

    # create synapse
    learning_rate = 0.1
    syn = PreEventSynapse(learning_rate, tau_trace=20.e-3, A_minus=-0.15)

    # integration time step
    dt = 0.001

    weights = []
    bp = []
    for r in rates:
        syn.reset()
        protocol = RandomProtocol(duration=duration, rate=r)
        #protocol.display('../../results/learning-rule/Poisson/protocolRandom_r'+str(r)[:2]+'.pdf')
        syn.pre.compute_trains(protocol.spiketimes_pre, dt, duration)
        syn.post.compute_trains(protocol.spiketimes_post, dt, duration)
        bp.append(np.sum(syn.pre.train['burst'])/np.sum(syn.pre.train['event']))
        print(r, bp[-1])

        for t in range(int(duration/dt)):
            syn.evolve(t, dt)
        #print(syn)
        weights.append(syn.weight/dt)
        #  division by dt is because event and burst trains are not divided
        #  by dt in the rule (see preeventsynapse.py)
    return rates, weights

def analytical(eta, duration, A_plus, A_minus, burst_threshold, tau_pre):
    # select parameters for pairing protocol
    rates = np.linspace(1., 20., 50)
    p_0 = -A_minus/A_plus
    a = np.exp(-burst_threshold*rates)
    event_rates = rates*a
    weights = duration*eta*event_rates**2*tau_pre*(1. - a - p_0)

    return rates, weights


def burst_fraction(rates, burst_threshold):
    return 1-np.exp(-burst_threshold*rates)

if __name__ == '__main__':
    r, w = simulate()
    r_anal, w_anal = analytical(0.1, 400, 1., -0.15, 0.016, 0.020)
    #print(burst_fraction(np.linspace(1., 30., 10), 0.016))
    plt.plot(r,w, label='simul.')
    plt.plot(r_anal, w_anal, label='anal.')
    plt.plot(r, np.zeros(r.shape), 'k--', lw=1)
    plt.legend()
    plt.tight_layout()
    plt.show()