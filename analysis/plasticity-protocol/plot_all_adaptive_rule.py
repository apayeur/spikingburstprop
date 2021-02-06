# python 3.6.2
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

import simulate_burstpoissonprotocol as bpp
import simulate_periodicprotocol as pp
import simulate_randomprotocol as rp

plt.style.use('../thesis_mplrc.dms')

PLOTPATH = os.path.join('..', '..', 'results', 'learning-rule')
if not os.path.exists(PLOTPATH):
    os.makedirs(PLOTPATH)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate all pairing protocols
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters shared by all protocols
burst_threshold = 0.016
tau_pre = 0.050
learning_rate = 0.1
alpha = 15.

# Burst-Poisson protocol
duration = 100.
bps = np.arange(0., 0.5, 0.05)
w_bpp, std_bpp = bpp.simulate(bps, learning_rate, duration, alpha, burst_threshold, tau_pre, 5.)
w_bpp_other_ER, std_bpp_other_ER = bpp.simulate(bps, learning_rate, duration, alpha, burst_threshold, tau_pre, 10.)

# Periodic protocol
freqs = np.linspace(1., 100., 20)  # frequency in Hz
f, w_pp = pp.simulate(freqs, learning_rate, alpha, burst_threshold, tau_pre)

# Poisson protocol
rates = np.linspace(1., 50., 10)
duration = 100.
w_rp, std_rp = rp.simulate(rates, learning_rate, duration, alpha, burst_threshold, tau_pre, (5., 0.30*5.))
w_rp2, std_rp2 = rp.simulate(rates, learning_rate, duration, alpha, burst_threshold, tau_pre, (5., 0.50*5.))


# ~~~~~~~~
# Plotting
# ~~~~~~~~
fig, axs = plt.subplots(2, 2, figsize=(88./25.4*4.25/3.75, 88./25.4))
n_SD = 1  # width of error bars for Poisson and burst-Poisson protocols

# Empty subplot occupied by the schematic
axs[0, 0].set_title('Learning rule')
sns.despine(ax=axs[0, 0], bottom=True, left=True)
axs[0, 0].set_xticks([], [])
axs[0, 0].set_yticks([], [])

# Periodic protocol
axs[0, 1].plot(freqs, w_pp, color='black')
axs[0, 1].plot(freqs, np.zeros(len(freqs)), 'k--', lw=0.5)
axs[0, 1].set_xticks([1, 100])
axs[0, 1].set_yticks([-0.3, 0, 0.3])
axs[0, 1].set_xlabel('Pairing frequency [Hz]')
axs[0, 1].set_ylabel('$\Delta W$')
axs[0, 1].set_title('Periodic protocol')
sns.despine(ax=axs[0, 1])

# Poisson protocol
axs[1, 0].plot(rates, w_rp, color='black', label=r'$\overline{P}(0)$ = 30%')
axs[1, 0].fill_between(rates, w_rp - n_SD*std_rp, w_rp + n_SD*std_rp, color='black', alpha=0.25, lw=0)
axs[1, 0].plot(rates, w_rp2, color='grey', label=r'$\overline{P}(0)$ = 50%')
axs[1, 0].fill_between(rates, w_rp2 - n_SD*std_rp2, w_rp2 + n_SD*std_rp2, color='grey', alpha=0.25, lw=0)
axs[1, 0].plot(rates, np.zeros(len(rates)), 'k--', lw=0.5)
axs[1, 0].set_xticks([1, 50])
axs[1, 0].set_xlabel('Rate [Hz]')
axs[1, 0].set_ylabel('$\Delta W$')
axs[1, 0].legend(loc='center left', frameon=False, framealpha=0.5, ncol=1, handlelength=2, borderaxespad=0.1)
axs[1, 0].set_title('Poisson protocol')
sns.despine(ax=axs[1, 0])

# Burst-Poisson protocol
axs[1, 1].plot(100*bps, w_bpp, color='black', label="ER = 5 Hz")
axs[1, 1].fill_between(100*bps, w_bpp - n_SD*std_bpp, w_bpp + n_SD*std_bpp, color='black', alpha=0.25, lw=0)
axs[1, 1].plot(100*bps, w_bpp_other_ER, color='grey', label="ER = 10 Hz")
axs[1, 1].fill_between(100*bps, w_bpp_other_ER - n_SD*std_bpp_other_ER,
                       w_bpp_other_ER + n_SD*std_bpp_other_ER, color='grey', alpha=0.25, lw=0)
axs[1, 1].plot(100*bps, np.zeros(len(bps)), 'k--', lw=0.5)
axs[1, 1].set_xticks([0, 50])
axs[1, 1].set_yticks([0, 2])
axs[1, 1].set_xlabel('Burst probability [%]')
axs[1, 1].set_ylabel('$\Delta W$')
axs[1, 1].legend(loc='lower right')
axs[1, 1].set_title('Burst-Poisson protocol')
sns.despine(ax=axs[1, 1])

fig.align_ylabels(axs[:, 1])  # align y labels of right-hand panels

plt.tight_layout()
plt.savefig(PLOTPATH + '/WeightALLProtocols_AdaptiveRule.pdf')
plt.close()
