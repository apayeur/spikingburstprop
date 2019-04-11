import matplotlib.pyplot as plt
from custom_colors import custom_colors
import seaborn as sns
from pairingprotocols import PoissonBurstProtocol
from preeventsynapse import PreEventSynapse
import numpy as np

plt.style.use('/Users/alexandrepayeur/Documents/Recherches/Burstprop/Codes/thesis_mplrc.dms')


"""
Description:
The plasticity protocol consists in pre- and post-synaptic Poisson event-burst
processes. Namely, a neuron fires an event, and a burst is then evoked with probability p.
"""


# select parameters for pairing protocol
burst_probs = np.arange(0., 0.5, 0.05)
duration = 500.  # in seconds

# create synapse
learning_rate = 0.1
syn = PreEventSynapse(learning_rate, tau_trace=20.e-3)

# integration time step
dt = 0.001

weights = []
for p in burst_probs:
    syn.reset()
    protocol = PoissonBurstProtocol(p=p, duration=duration, rate_event=5.)
    protocol.display('../../results/learning-rule/Poisson/protocol_p'+str(p)[:4]+'.pdf')
    syn.pre.compute_trains(protocol.spiketimes_pre, dt, duration)
    syn.post.compute_trains(protocol.spiketimes_post, dt, duration)

    for t in range(int(duration/dt)):
        syn.evolve(t, dt)
    print(syn)
    weights.append(syn.weight/dt)
    #  division by dt is because event and burst trains are not divided
    #  by dt in the rule (see preeventsynapse.py)


# plotting
plt.figure(figsize=(1.5*1.6, 1.5))
plt.plot(100*burst_probs, weights, color='black')
plt.plot(100*burst_probs, np.zeros(len(burst_probs)), 'k--', lw=1)
sns.despine()
plt.xlabel('Postsyn. BP [\%]')
plt.ylabel('$\Delta W$')
plt.tight_layout()
plt.savefig('../../results/learning-rule/Poisson/Weight_vs_BurstProb.pdf')
plt.close()
