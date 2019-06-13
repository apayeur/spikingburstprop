import matplotlib.pyplot as plt
import numpy as np
import sys
import seaborn as sns

wsum = np.loadtxt('../../../data/learning-rule/EBCP/tau_pre20.e-3/d00.2/w00.045/plasticity_rule.0.wsum')


plt.figure(figsize=(4,4/1.6))
plt.plot(wsum[:,0], wsum[:,1], color='black', lw=3)
plt.xlabel('Time [s]')
plt.ylabel('summed weights')
plt.tight_layout()
sns.despine()
plt.savefig('Wsum.pdf')
plt.close()
