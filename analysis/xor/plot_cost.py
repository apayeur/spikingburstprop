import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../')
plt.style.use('../thesis_mplrc.dms')
from utils_plotting import remove_xticklabel
from utils_plotting import custom_colors


###     DEBUG    ###
cost = np.loadtxt('../../data/xor/very_good_xor/cost2.dat')
epoch_cost = np.zeros((cost.shape[0]//4, 3))

for i in range(cost.shape[0]//4):
    epoch_cost[i,0] = i
    epoch_cost[i, 1] = np.sum(cost[4*i:4*(i+1),1])
    epoch_cost[i, 2] = np.sum(cost[4*i:4*(i+1),2])


# Cost for all examples consecutively
plt.figure(figsize=(3,3/1.5))
plt.plot(cost[:,0], cost[:,1], 'k', lw=0.2)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Example \#')
plt.ylabel('Cost')
plt.tight_layout()
plt.savefig('../../results/xor/Cost.pdf')
plt.close()


# Cost vs epoch
plt.figure(figsize=(3, 4/1.5))
plt.plot(epoch_cost[:,0], epoch_cost[:,1], 'k')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Epoch')
plt.ylabel('Cost')
plt.tight_layout()
plt.savefig('../../results/xor/Epoch_Cost.pdf')
plt.close()


# Cost for each type of inputs individually
fig = plt.figure(2, figsize=(4, 4/1.5))
ax1 = fig.add_subplot(221)
ax1.plot(epoch_cost[:,0], cost[0::4, 1])
ax1.set_ylabel('Cost for (0, 0)')
remove_xticklabel(ax1)

ax2 = fig.add_subplot(222)
ax2.plot(epoch_cost[:,0], cost[1::4, 1])
ax2.set_ylabel('Cost for (1, 0)')
remove_xticklabel(ax2)

ax3 = fig.add_subplot(223)
ax3.plot(epoch_cost[:,0], cost[2::4, 1])
ax3.set_xlabel('Epoch')
ax3.set_ylabel('Cost for (0, 1)')

ax4 = fig.add_subplot(224)
ax4.plot(epoch_cost[:,0], cost[3::4, 1])
ax4.set_ylabel('Cost for (1, 1)')
ax4.set_xlabel('Epoch')

plt.tight_layout()
plt.savefig('../../results/xor/Cost_EachExample.pdf')
plt.close()
