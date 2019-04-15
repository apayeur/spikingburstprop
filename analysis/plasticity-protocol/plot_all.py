# python 3.6.2
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

import simulate_burstpoissonprotocol as bpp
import simulate_periodicprotocol as pp
import simulate_randomprotocol as rp


plt.style.use('../thesis_mplrc.dms')


# simulate all pairing protocols
(burst_probs, w_bpp) = bpp.simulate()
(freqs, w_pp) = pp.simulate()
(rates, w_rp) = rp.simulate()

# analytical results
freqs_anal, w_pp_anal = pp.analytical(0.1, 1., -0.15, 0.016)
(burst_probs_anal, w_bpp_anal) = bpp.analytical(0.1, 400, 5., 1, -0.15, 0.016, 0.020, 0.020)
(rates_anal, w_rp_anal) = rp.analytical(0.1, 400, 1., -0.15, 0.016, 0.020)


plt.figure(figsize=(1.5*1.6, 1.5*3))

plt.subplot(311)
plt.plot(freqs, w_pp, color='black')
plt.plot(freqs_anal, w_pp_anal, color='grey', lw=1)
plt.plot(freqs, np.zeros(len(freqs)), 'k--', lw=0.5)
sns.despine()
plt.xlabel('Pairing frequency [Hz]')
plt.ylabel('$\Delta W$')
plt.title('Periodic protocol')

plt.subplot(312)
plt.plot(rates, w_rp, color='black')
plt.plot(rates_anal, w_rp_anal, color='grey', lw=1)
plt.plot(rates, np.zeros(len(rates)), 'k--', lw=0.5)
sns.despine()
plt.xlabel('Rate [Hz]')
plt.ylabel('$\Delta W$')
plt.title('Poisson protocol')

plt.subplot(313)
plt.plot(100*burst_probs, w_bpp, color='black')
plt.plot(100*burst_probs_anal, w_bpp_anal, color='grey', lw=1)
plt.plot(100*burst_probs, np.zeros(len(burst_probs)), 'k--', lw=0.5)
sns.despine()
plt.xlabel('Postsyn. BP [\%]')
plt.ylabel('$\Delta W$')
plt.title('Burst-Poisson protocol')

plt.tight_layout()
plt.savefig('../../results/learning-rule/WeightALLProtocols.pdf')
plt.close()
