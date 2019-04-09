import matplotlib.pyplot as plt
from custom_colors import custom_colors
import seaborn as sns
from pairingprotocols import PeriodicPairing
from preeventsynapse import PreEventSynapse
import numpy as np

plt.style.use('/Users/alexandrepayeur/Documents/Recherches/Burstprop/Codes/thesis_mplrc.dms')


"""
Description:
The plasticity protocol consists in 15 bursts of 5 pairings, assuming a long
time in between bursts. This way, we can simply consider 5 pairings
and multiply the weight change by 15 at the end (since there is no weight decay).
"""


# select parameters for pairing protocol
freqs = np.arange(2, 80., 2.)  # frequency in Hz
Delta = 0.010  # t_post - t_pre
number_of_pairings = 5

# create synapse
syn = PreEventSynapse(0.1, burst_def=0.025)

# integration time step
dt = 0.001

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

# plotting
plt.figure(figsize=(2.*1.6, 2.))
plt.plot(freqs, weights, color='black')
plt.plot(freqs, np.zeros(len(freqs)), 'k--', lw=1)
sns.despine()
plt.xlabel('Pairing frequency [Hz]')
plt.ylabel('$\Delta W$')
plt.title('$\Delta t = $' + str(Delta*1000)[:2] + 'ms')
plt.tight_layout()
plt.savefig('../../results/learning-rule/Periodic/Weight_vs_Freq_Delta'+str(Delta*1000)[:2]+'.pdf')
plt.close()


