# python 3.6.2
import matplotlib.pyplot as plt
#import seaborn as sns
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
duration = 100.
bps = np.arange(0., 0.5, 0.05)
alpha = 30.
w_bpp = bpp.simulate(bps, learning_rate, duration, alpha, burst_threshold, tau_pre, 5.)
w_bpp_other_ER = bpp.simulate(bps, learning_rate, duration, alpha, burst_threshold, tau_pre, 10.)

# Periodoc protocol
freqs = np.arange(2, 80., 2.)  # frequency in Hz
alpha = 10.
f, w_pp = pp.simulate(freqs, alpha)

# Poisson protocol
rates = np.linspace(1, 50, 10)
duration = 60.
r, w_rp, _, _ = rp.simulate(rates, learning_rate, duration, 30., burst_threshold, tau_pre, (5., 0.30*5.))
r, w_rp2, _, _ = rp.simulate(rates, learning_rate, duration, 30., burst_threshold, tau_pre, (5., 0.50*5.))


# --- plotting --- #
plt.figure(figsize=(4.25, 3.75))

plt.subplot(222)
plt.title('Learning rule')


plt.subplot(222)
plt.plot(freqs, w_pp, color='black')
plt.plot(freqs, np.zeros(len(freqs)), 'k--', lw=0.5)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Pairing frequency [Hz]')
plt.ylabel('$\Delta W$')
plt.title('Periodic protocol')

plt.subplot(223)
plt.plot(rates, w_rp, color='black', label=r'$\bar{P}(0) = 30\%$')
plt.plot(rates, w_rp2, color='grey', label=r'$\bar{P}(0) = 50\%$')
plt.plot(rates, np.zeros(len(rates)), 'k--', lw=0.5)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Rate [Hz]')
plt.ylabel('$\Delta W$')
plt.legend(loc = 'upper left')
plt.title('Poisson protocol')

plt.subplot(224)
plt.plot(100*bps, w_bpp, color='black', label="ER = 5Hz")
plt.plot(100*bps, w_bpp_other_ER, color='grey', label="ER = 10Hz")
plt.plot(100*bps, np.zeros(len(bps)), 'k--', lw=0.5)
plt.yticks([-4,0,4])
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Burst probability [\%]')
plt.ylabel('$\Delta W$')
plt.legend(loc = 'lower right')
plt.title('Burst-Poisson protocol')

plt.tight_layout()
plt.savefig('../../results/learning-rule/WeightALLProtocols_AdaptiveRule.pdf')
plt.close()
