import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../')
plt.style.use('../thesis_mplrc.dms')
from utils_plotting import remove_xticklabel
from utils_plotting import custom_colors

cost = np.loadtxt('../../data/xor/cost2.dat')
epoch_cost = np.zeros((cost.shape[0]//4, 3))

for i in range(cost.shape[0]//4):
    epoch_cost[i,0] = i
    epoch_cost[i, 1] = np.sum(cost[4*i:4*(i+1),1])
    epoch_cost[i, 2] = np.sum(cost[4*i:4*(i+1),2])


plt.figure(figsize=(3,3/1.5))
plt.plot(cost[:,0], cost[:,1], 'k', lw=0.2)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Example \#')
plt.ylabel('Cost')
plt.tight_layout()
plt.savefig('../../results/xor/Cost.pdf')
plt.close()


plt.figure(figsize=(3,3/1.5))
plt.plot(epoch_cost[:,0], epoch_cost[:,1], 'k')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Epoch')
plt.ylabel('Cost')
plt.tight_layout()
plt.savefig('../../results/xor/Epoch_Cost.pdf')
plt.close()
