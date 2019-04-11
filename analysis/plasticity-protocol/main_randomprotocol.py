import matplotlib.pyplot as plt
from custom_colors import custom_colors
import seaborn as sns
from pairingprotocols import RandomProtocol
from preeventsynapse import PreEventSynapse
import numpy as np

plt.style.use('/Users/alexandrepayeur/Documents/Recherches/Burstprop/Codes/thesis_mplrc.dms')


"""
Description:
The plasticity protocol consists in pre- and post-synaptic Poisson processes.
"""


# select parameters for pairing protocol
rates = np.linspace(1., 30., 10)

# create synapse
learning_rate = 0.1
syn = PreEventSynapse(learning_rate, tau_trace=20.e-3, A_minus=-0.15)

# integration time step
dt = 0.001

weights = []
for r in rates:
    syn.reset()
    protocol = RandomProtocol(spike_number=2000, rate=r)
    protocol.display('../../results/learning-rule/Poisson/protocolRandom_r'+str(r)[:2]+'.pdf')
    duration = max(np.max(protocol.spiketimes_pre), np.max(protocol.spiketimes_post)) + 0.1
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
plt.plot(rates, weights, color='black')
plt.plot(rates, np.zeros(len(rates)), 'k--', lw=1)
sns.despine()
plt.xlabel('Rate [Hz]')
plt.ylabel('$\Delta W$')
plt.tight_layout()
plt.savefig('../../results/learning-rule/Poisson/Weight_vs_Rates_RandomProtocol.pdf')
plt.close()
