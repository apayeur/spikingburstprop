import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../../analysis')
plt.style.use('../../../analysis/thesis_mplrc.dms')
from utils_plotting import remove_xticklabel
from utils_plotting import custom_colors
plt.rcParams["axes.labelsize"] = 8
plt.rcParams["legend.fontsize"] = 6
plt.rcParams["legend.borderpad"] = 0.25
plt.rcParams["legend.labelspacing"] = 0.1
plt.rcParams["xtick.labelsize"] = 7
plt.rcParams["ytick.labelsize"] = 7



# data folder
DATAPATH = '../data/'

# import data
lb = 0
ub = int(50./0.1)
#[lb:ub, :]
p2 = np.loadtxt(DATAPATH + 'p_output.dat')[lb:ub, :]
p2_bar = np.loadtxt(DATAPATH + 'p_bar_output.dat')[lb:ub, :]
p1 = np.loadtxt(DATAPATH + 'p_hidden.dat')[lb:ub, :]
p1_bar = np.loadtxt(DATAPATH + 'p_bar_hidden.dat')[lb:ub, :]

e2 = np.loadtxt(DATAPATH + 'e_output.dat')[lb:ub, :]
e1 = np.loadtxt(DATAPATH + 'e_hidden.dat')[lb:ub, :]
e0 = np.loadtxt(DATAPATH + 'e_input.dat')[lb:ub, :]
d = np.loadtxt(DATAPATH + 'teacher.dat')[lb:ub, :]

W = np.loadtxt(DATAPATH + 'example_weights.dat')[lb:ub, :]

# time
t = 25*64 + np.arange(0, 50., 0.1)
xtic = [t[0], t[0]+25, t[0]+50]

# plotting
plt.figure(1, figsize=(4, 3.5))

# output layer
ax1 = plt.subplot(521)
ax1.plot(t, e2[:, 5], color=custom_colors['blue'], lw=1)
ax1.set_ylim([0,0.2])
ax1.set_yticks([0,0.2])
ax1.set_ylabel(r'$e_5^{(2)}$', color=custom_colors['blue'])
ax1.spines['top'].set_visible(False)
ax1.spines['bottom'].set_color(custom_colors['pink'])
ax1.spines['left'].set_color(custom_colors['pink'])
ax1.set_xticks(xtic)
remove_xticklabel(ax1)

ax1_tw = ax1.twinx()
ax1_tw.plot(t, d[:, 5], '--', color=custom_colors['bluish_green'], lw=0.5)
ax1_tw.set_ylabel('Teacher$_5$', color=custom_colors['bluish_green'])
ax1_tw.spines['top'].set_visible(False)
remove_xticklabel(ax1_tw)

ax4 = plt.subplot(523)
ax4.plot(t, W[:, 1], color='black', lw=1)
ax4.set_ylabel('$W_{5,100}^{(2)}$')
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.legend(loc='best')
ax4.set_xticks(xtic)
#ax4.set_yticks([-0.0788, -0.0793])
remove_xticklabel(ax4)

ax5 = plt.subplot(525)
ax5.plot(t, e1[:, 100], color=custom_colors['blue'], lw=1)
ax5.set_ylabel(r'$e_{100}^{(1)}$', color=custom_colors['blue'])
ax5.spines['top'].set_visible(False)
ax5.spines['right'].set_visible(False)
ax5.set_xticks(xtic)
remove_xticklabel(ax5)

ax8 = plt.subplot(527)
ax8.plot(t, W[:, 0], color='black', lw=1)
ax8.set_ylabel('$W_{100,300}^{(1)}$')
#ax8.set_ylim([0.00263,0.00265])
ax8.spines['top'].set_visible(False)
ax8.spines['right'].set_visible(False)
ax8.legend(loc='best')
ax8.set_xticks(xtic)
remove_xticklabel(ax8)

# input layer
ax9 = plt.subplot(529)
ax9.plot(t, e0[:, 300], color=custom_colors['blue'], lw=1)
ax9.set_ylabel(r'$e_{300}^{(0)}$', color=custom_colors['blue'])
ax9.set_xlabel('Time [s]')
ax9.spines['top'].set_visible(False)
ax9.spines['right'].set_visible(False)
ax9.set_xticks(xtic)

# BP
ax2 = plt.subplot(522)
ax2.plot(t, p2[:, 5], color=custom_colors['red'], lw=0.5, label='$p_5^{(2)}$')
ax2.plot(t, p2_bar[:, 5], '--',  color=custom_colors['red'], lw=0.5, label=r'$\bar{p}_5^{(2)}$')
ax2.set_ylabel('Burst prob.', color=custom_colors['red'])
ax2.set_yticks([0,1])
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.legend(loc='upper right', bbox_to_anchor=(1, 1.4))
ax2.set_xticks(xtic)
remove_xticklabel(ax2)

ax3 = plt.subplot(524)
ax3.plot(t, p2[:, 5]-p2_bar[:, 5], color=custom_colors['red'], lw=0.5, label='$p_5^{(2)}$')
ax3.plot(t, np.zeros(p2[:, 5].shape), '--', lw=0.2, color='black')
#ax3.plot(t, p2_bar[:, 5], '--', color=custom_colors['red'], lw=1, label=r'$\bar{p}_5^{(2)}$')
ax3.text(0.5, 1.2, r'$\delta p_5^{(2)} = p_5^{(2)} - \bar{p}_5^{(2)}$', color=custom_colors['red'],
         horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, fontsize=6)
ax3.set_ylabel(r'$\delta p_5^{(2)}$', color=custom_colors['red'])
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.set_xticks(xtic)
remove_xticklabel(ax3)

ax6 = plt.subplot(526)
ax6.plot(t, p1[:, 100], color=custom_colors['red'], lw=0.5, label='$p_{100}^{(1)}$')
ax6.plot(t, p1_bar[:, 100], '--',  color=custom_colors['red'], lw=0.5, label=r'$\bar{p}_{100}^{(1)}$')
ax6.set_ylabel('Burst prob.', color=custom_colors['red'])
#ax6.set_yticks([0,1])
ax6.spines['top'].set_visible(False)
ax6.spines['right'].set_visible(False)
ax6.legend(loc='best')
ax6.set_xticks(xtic)
remove_xticklabel(ax6)

ax7 = plt.subplot(528)
ax7.plot(t, p1[:, 100]-p1_bar[:, 100], color=custom_colors['red'], lw=0.5, label='$p_{100}^{(1)}$')
ax7.plot(t, np.zeros(p1[:, 100].shape), '--', lw=0.2, color='black')
ax7.set_ylabel(r'$\delta p_{100}^{(1)}$', color=custom_colors['red'])
ax7.text(0.5, 1.2, r'$\delta p_{100}^{(1)} = p_{100}^{(1)} - \bar{p}_{100}^{(1)}$', color=custom_colors['red'],
         horizontalalignment='center', verticalalignment='center', transform=ax7.transAxes, fontsize=6)
ax7.spines['top'].set_visible(False)
ax7.spines['right'].set_visible(False)
ax7.set_xticks(xtic)
remove_xticklabel(ax7)

plt.tight_layout()
plt.savefig('../results/RateModel.pdf')
plt.close()



