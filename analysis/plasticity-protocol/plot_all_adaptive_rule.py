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
tau_pre = 0.050
learning_rate = 0.1
alpha = 15.


# Burst-Poisson protocol
duration = 100.
bps = np.arange(0., 0.5, 0.05)
w_bpp = bpp.simulate(bps, learning_rate, duration, alpha, burst_threshold, tau_pre, 5.)
w_bpp_other_ER = bpp.simulate(bps, learning_rate, duration, alpha, burst_threshold, tau_pre, 10.)

# Periodic protocol
freqs = np.linspace(1., 100., 20)  # frequency in Hz
f, w_pp = pp.simulate(freqs, learning_rate, alpha, burst_threshold, tau_pre)

# Poisson protocol
rates = np.linspace(1., 50., 10)
duration = 60.
r, w_rp, _, _ = rp.simulate(rates, learning_rate, duration, alpha, burst_threshold, tau_pre, (5., 0.30*5.))
r, w_rp2, _, _ = rp.simulate(rates, learning_rate, duration, alpha, burst_threshold, tau_pre, (5., 0.50*5.))


# --- plotting --- #
fig, axs = plt.subplots(2, 2, figsize=(88./25.4*4.25/3.75, 88./25.4))

axs[0, 0].set_title('Learning rule')
sns.despine(ax=axs[0, 0], bottom=True, left=True)
axs[0, 0].set_xticks([], [])
axs[0, 0].set_yticks([], [])

axs[0, 1].plot(freqs, w_pp, color='black')
axs[0, 1].plot(freqs, np.zeros(len(freqs)), 'k--', lw=0.5)
axs[0, 1].set_xticks([1, 100])
axs[0, 1].set_yticks([-0.3, 0, 0.3])
axs[0, 1].set_xlabel('Pairing frequency [Hz]')
axs[0, 1].set_ylabel('$\Delta W$')
axs[0, 1].set_title('Periodic protocol')
sns.despine(ax=axs[0, 1])

axs[1, 0].plot(rates, w_rp, color='black', label=r'$\overline{P}(0)$ = 30%')
axs[1, 0].plot(rates, w_rp2, color='grey', label=r'$\overline{P}(0)$ = 50%')
axs[1, 0].plot(rates, np.zeros(len(rates)), 'k--', lw=0.5)
axs[1, 0].set_xticks([1, 50])
axs[1, 0].set_xlabel('Rate [Hz]')
axs[1, 0].set_ylabel('$\Delta W$')
axs[1, 0].legend(loc='lower right', frameon=False, framealpha=0.5, ncol=1, handlelength=1, borderaxespad=0.1, prop={'size': 7})
axs[1, 0].set_title('Poisson protocol')
sns.despine(ax=axs[1, 0])

axs[1, 1].plot(100*bps, w_bpp, color='black', label="ER = 5 Hz")
axs[1, 1].plot(100*bps, w_bpp_other_ER, color='grey', label="ER = 10 Hz")
axs[1, 1].plot(100*bps, np.zeros(len(bps)), 'k--', lw=0.5)
axs[1, 1].set_xticks([0, 50])
axs[1, 1].set_yticks([0, 2])
axs[1, 1].set_xlabel('Burst probability [%]')
axs[1, 1].set_ylabel('$\Delta W$')
axs[1, 1].legend(loc='lower right')
axs[1, 1].set_title('Burst-Poisson protocol')
sns.despine(ax=axs[1, 1])

# align y labels of right-hand panels
fig.align_ylabels(axs[:, 1])

plt.tight_layout()
plt.savefig('../../results/learning-rule/WeightALLProtocols_AdaptiveRule.pdf')
plt.close()
