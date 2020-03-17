# python 3.6.2
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

import simulate_randomprotocol_adaptive_rule as rp

plt.style.use('../thesis_mplrc.dms')


burst_threshold = 0.016
tau_pre = 0.050
learning_rate = 0.1
alpha = 15.
duration = 60.
rates = np.linspace(1., 50., 50)

'''
def find_zero(x):
    position_of_min = np.argmin(np.abs(x))
    if DeltaW[position_of_min] > 0:
        slope =
        return DeltaW[position_of_min]

pbar0s = np.linspace(0, 50, 25)
threshold = []
for p in pbar0s:
    r, w_rp, _, _ = rp.simulate(rates, learning_rate, duration, alpha, burst_threshold, tau_pre, (5., p*5.))
    threshold.append(find_zero(w_rp))

plt.figure(figsize=(4, 3))
plt.plot(r, threshold)
plt.ylabel('Threshold')
plt.xlabel('Rate [Hz]')
sns.despine(trim=True)
plt.show()
'''


def DeltaW(r, p0, T):
    E = r*np.exp(-r*0.016)
    return alpha*learning_rate*tau_pre*(1 - np.exp(-r*0.016) - p0)*(1 - np.exp(-T/alpha))*E*E


plt.figure(figsize=(3, 3/1.6))
plt.plot(rates, DeltaW(rates, 0.3, 60)/rates, color='black', label='$\overline{P}(0)$ = 30%')
plt.plot(rates, DeltaW(rates, 0.4, 60)/rates, '--', color='black', label='$\overline{P}(0)$ = 40%')
plt.plot(rates, DeltaW(rates, 0.5, 60)/rates, color='grey', label='$\overline{P}(0)$ = 50%')
plt.plot(rates, np.zeros(rates.shape), '--', color='black', lw=1)
plt.xlabel('Rate [Hz]')
plt.ylabel('$\phi$')
sns.despine(trim=True)
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig('../../results/learning-rule/PhiBCM.png')
plt.close()


plt.figure(figsize=(3, 3/1.6))
p0s = np.linspace(0,0.5,50)
plt.plot(rates, (1 - np.exp(-0.016*rates)), color='black')
plt.ylabel('Threshold')
plt.xlabel('Rate [Hz]')
sns.despine(trim=True)
plt.tight_layout()
plt.savefig('../../results/learning-rule/ThresholdTheory.png')
plt.close()
