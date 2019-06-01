# python 3.6.2
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

import simulate_burstpoissonprotocol_adaptive_rule as bpp
import simulate_periodicprotocol_adaptive_rule as pp
import simulate_randomprotocol_adaptive_rule as rp

plt.style.use('../thesis_mplrc.dms')

# --- simulate all pairing protocols --- #
# Parameters shared by all protocols
burst_threshold = 0.016
tau_pre = 0.020
learning_rate = 0.1

# Burst-Poisson protocol
duration = 60.
bps = np.arange(0., 0.5, 0.05)
alpha = 30.
w_bpp = bpp.simulate(bps, learning_rate, duration, alpha, burst_threshold, tau_pre)

# Periodoc protocol
freqs = np.arange(2, 80., 2.)  # frequency in Hz
alpha = 10.
f, w_pp = pp.simulate(freqs, alpha)

# Poisson protocol
rates = np.linspace(1, 50, 10)
duration = 60.
r, w_rp, _, _ = rp.simulate(rates, learning_rate, duration, 30., burst_threshold, tau_pre)


# --- plotting --- #
plt.figure(figsize=(1.5*1.6, 1.5*3))

plt.subplot(311)
plt.plot(freqs, w_pp, color='black')
plt.plot(freqs, np.zeros(len(freqs)), 'k--', lw=0.5)
sns.despine()
plt.xlabel('Pairing frequency [Hz]')
plt.ylabel('$\Delta W$')
plt.title('Periodic protocol')

plt.subplot(312)
plt.plot(rates, w_rp, color='black')
plt.plot(rates, np.zeros(len(rates)), 'k--', lw=0.5)
sns.despine()
plt.xlabel('Rate [Hz]')
plt.ylabel('$\Delta W$')
plt.title('Poisson protocol')

plt.subplot(313)
plt.plot(100*bps, w_bpp, color='black')
plt.plot(100*bps, np.zeros(len(bps)), 'k--', lw=0.5)
sns.despine()
plt.xlabel('Initial BP [\%]')
plt.ylabel('$\Delta W$')
plt.title('Burst-Poisson protocol')

plt.tight_layout()
plt.savefig('../../results/learning-rule/WeightALLProtocols_AdaptiveRule.pdf')
plt.close()
