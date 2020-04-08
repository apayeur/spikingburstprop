import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import sys
sys.path.append('../')
plt.style.use('../thesis_mplrc.dms')
plt.rcParams["axes.labelsize"] = 9
from utils_plotting import remove_xticklabel
from utils_plotting import custom_colors


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-durex", type=float, help="duration of example")
parser.add_argument("-alpha", type=float, help="moving average time scale")
parser.add_argument("-numex", type=int, help="number of training examples")
args = parser.parse_args()


# import data
data_file_prefix = '../../data/xor/xor.0.brate_'
datadir = '../../data/xor/'

output_rates = np.loadtxt(data_file_prefix + 'output')
hidden1_rates = np.loadtxt(data_file_prefix + 'hidden1')
hidden2_rates = np.loadtxt(data_file_prefix + 'hidden2')
input1_rates = np.loadtxt(data_file_prefix + 'input1')
input2_rates = np.loadtxt(data_file_prefix + 'input2')
cost_epoch = np.loadtxt(datadir + 'cost_epoch.dat')

binSize = output_rates[1,0] - output_rates[0,0]
print(binSize)
after_learning = [output_rates.shape[0]-int(4*args.durex/binSize)+1, output_rates.shape[0]+1]
before_learning = [int(3*args.alpha/binSize)+1, int(3*args.alpha/binSize)+int(4*args.durex/binSize)+1]
during_learning = [int(3*args.alpha/binSize)+int(args.numex//2*args.durex/binSize)+1,
                   int(3*args.alpha/binSize)+int(args.numex//2*args.durex/binSize)+int(4*args.durex/binSize)+1]


er_trace = np.loadtxt(datadir + 'xor0.0.trevent')
br_trace = np.loadtxt(datadir + 'xor0.0.trburst')
for i in range(1,50):
    er_tmp = np.loadtxt(datadir + 'xor'+str(i)+'.0.trevent')
    br_tmp = np.loadtxt(datadir + 'xor'+str(i)+'.0.trburst')
    er_trace[:,1] += er_tmp[:,1]
    br_trace[:,1] += br_tmp[:,1]
er_trace[:,1] /= 50.
br_trace[:,1] /= 50.
er_trace = er_trace
br_trace = br_trace

data_file_prefix = '../../data/xor/xor.0.wsum_'
wsum_in1_to_hid1 = np.loadtxt(data_file_prefix + 'in1_to_hid1')
wsum_in1_to_hid2 = np.loadtxt(data_file_prefix + 'in1_to_hid2')
wsum_in2_to_hid1 = np.loadtxt(data_file_prefix + 'in2_to_hid1')
wsum_in2_to_hid2 = np.loadtxt(data_file_prefix + 'in2_to_hid2')
wsum_hid1_to_out = np.loadtxt(data_file_prefix + 'hid1_to_out')
wsum_hid2_to_out = np.loadtxt(data_file_prefix + 'hid2_to_out')

e_min = 2.
e_max = 10.
e_th = e_min + 0.5*(e_max - e_min)


def plot_section(lim_seg, moment):

    # select data segments
    r_out = output_rates[lim_seg[0]:lim_seg[1]+1, :]
    r_hid1 = hidden1_rates[lim_seg[0]:lim_seg[1]+1, :]
    r_hid2 = hidden2_rates[lim_seg[0]:lim_seg[1]+1, :]
    r_inp1 = input1_rates[lim_seg[0]:lim_seg[1]+1, :]
    r_inp2 = input2_rates[lim_seg[0]:lim_seg[1]+1, :]

    er_est = er_trace[lim_seg[0]:lim_seg[1]+1, :]
    br_est = br_trace[lim_seg[0]:lim_seg[1]+1, :]

    xtik = list((lim_seg[0]-1)*binSize + args.durex*np.arange(5))
    print(xtik)

    # plotting rates and burst probabilities
    fig = plt.figure(1, figsize=(5.5, 6/1.5))
    gs1 = gridspec.GridSpec(3, 2)

    #  output layer
    ax1 = fig.add_subplot(gs1[0, 1])
    ax1.plot(r_out[:, 0], r_out[:, 2], color=custom_colors['blue'], lw=1, label='ER')
    ax1.plot(r_out[:, 0], r_out[:, 1], color=custom_colors['orange'], lw=1, label='BR')
    ax1.plot(r_out[:, 0], e_min*np.ones((r_out.shape)), 'k--', lw=0.5)
    ax1.plot(r_out[:, 0], e_th*np.ones((r_out.shape)), 'k--', lw=0.5)
    ax1.plot(r_out[:, 0], e_max*np.ones((r_out.shape)), 'k--', lw=0.5)
    ax1.set_ylabel('Output rates [Hz]')
    ax1.set_ylim([0,10])
    ax1.set_xticks(xtik)
    ax1.spines['top'].set_visible(False)
    ax1.legend(loc='best')

    ax1_tw = ax1.twinx()
    ax1_tw.plot(r_out[:, 0], 100*r_out[:, 1]/r_out[:, 2], color=custom_colors['red'], lw=1)
    ax1_tw.plot(er_est[:, 0], 100*br_est[:, 1]/er_est[:, 1], '--', color=custom_colors['red'], lw=1)
    ax1_tw.tick_params('y', labelcolor=custom_colors['red'])
    ax1_tw.set_yticks([0, 30, 60, 90])
    ax1_tw.set_ylabel('Postsyn. BP [\%]', color=custom_colors['red'])
    remove_xticklabel(ax1)


    #  hidden layer 1
    ax2 = fig.add_subplot(gs1[1, 0])
    ax2.plot(r_hid1[:, 0], r_hid1[:, 2], color=custom_colors['blue'], lw=1, label='ER')
    ax2.plot(r_hid1[:, 0], r_hid1[:, 1], color=custom_colors['orange'], lw=1, label='BR')
    ax2.plot(r_hid1[:, 0], e_min*np.ones((r_hid1.shape)), 'k--', lw=0.5)
    ax2.plot(r_hid1[:, 0], e_th*np.ones((r_hid1.shape)), 'k--', lw=0.5)
    ax2.plot(r_hid1[:, 0], e_max*np.ones((r_hid1.shape)), 'k--', lw=0.5)
    ax2.set_ylabel('Hidden1 rates [Hz]')
    ax2.set_xticks(xtik)
    ax2.set_ylim([0,10])
    ax2.spines['top'].set_visible(False)

    ax2_tw = ax2.twinx()
    ax2_tw.plot(r_hid1[:, 0], 100*r_hid1[:, 1]/r_hid1[:, 2], color=custom_colors['red'], lw=1)
    ax2_tw.tick_params('y', labelcolor=custom_colors['red'])
    ax2_tw.set_yticks([0, 30, 60, 90])
    ax2_tw.set_ylabel('Postsyn. BP [\%]', color=custom_colors['red'])
    remove_xticklabel(ax2)

    #  hidden layer 2
    ax3 = fig.add_subplot(gs1[1, 1])
    ax3.plot(r_hid2[:, 0], r_hid2[:, 2], color=custom_colors['blue'], lw=1, label='ER')
    ax3.plot(r_hid2[:, 0], r_hid2[:, 1], color=custom_colors['orange'], lw=1, label='BR')
    ax3.plot(r_hid2[:, 0], e_min*np.ones((r_hid2.shape)), 'k--', lw=0.5)
    ax3.plot(r_hid2[:, 0], e_th*np.ones((r_hid2.shape)), 'k--', lw=0.5)
    ax3.plot(r_hid2[:, 0], e_max*np.ones((r_hid2.shape)), 'k--', lw=0.5)
    ax3.set_ylabel('Hidden2 rates [Hz]')
    ax3.set_ylim([0,10])
    ax3.set_xticks(xtik)
    ax3.spines['top'].set_visible(False)

    ax3_tw = ax3.twinx()
    ax3_tw.plot(r_hid2[:, 0], 100*r_hid2[:, 1]/r_hid2[:, 2], color=custom_colors['red'], lw=1)
    ax3_tw.tick_params('y', labelcolor=custom_colors['red'])
    ax3_tw.set_yticks([0, 30, 60, 90])
    #ax3_tw.set_xticks((xtik))
    ax3_tw.set_ylabel('Postsyn. BP [\%]', color=custom_colors['red'])
    remove_xticklabel(ax3)


    #  input layer 1
    ax4 = fig.add_subplot(gs1[2, 0])
    ax4.plot(r_inp1[:, 0], r_inp1[:, 2], color=custom_colors['blue'], lw=1, label='ER')
    ax4.plot(r_inp1[:, 0], e_min*np.ones((r_inp1.shape)), 'k--', lw=0.5)
    ax4.plot(r_inp1[:, 0], e_max*np.ones((r_inp1.shape)), 'k--', lw=0.5)
    ax4.set_xticks(xtik)
    ax4.set_xlabel("Time [s]")
    ax4.set_ylabel('ER, input1 [Hz]')
    ax4.set_ylim([0,10])

    #  input layer 2
    ax5 = fig.add_subplot(gs1[2, 1])
    ax5.plot(r_inp2[:, 0], r_inp2[:, 2], color=custom_colors['blue'], lw=1, label='ER')
    #ax5.plot(r_hid2[:, 0], e_min*np.ones((r_hid2.shape)), 'k--', lw=0.5)
    #ax5.plot(r_hid2[:, 0], e_max*np.ones((r_hid2.shape)), 'k--', lw=0.5)
    ax5.set_xlabel("Time [s]")
    ax5.set_ylabel('ER, input2 [Hz]')
    ax5.set_ylim([0,10])
    ax5.set_xticks(xtik)

    plt.tight_layout()
    plt.savefig('../../results/xor/Activity_'+moment+'.pdf')
    plt.close()




plot_section(before_learning, 'Before')
plot_section(during_learning, 'DuringLearning')
plot_section(after_learning, 'After')

# plotting sum of weights
fig = plt.figure(2, figsize=(5, 5/1.5))
ax1 = fig.add_subplot(221)
ax1.plot(wsum_hid1_to_out[:,0], wsum_hid1_to_out[:,1])
ax1.set_ylabel('Average weight')
ax1.set_title('Hidden 1 to output')
#ax1.set_xticks([0]+tics)
remove_xticklabel(ax1)

ax2 = fig.add_subplot(222)
ax2.plot(wsum_hid2_to_out[:,0], wsum_hid2_to_out[:,1])
ax2.set_title('Hidden 2 to output')
#ax2.set_xticks([0]+tics)
remove_xticklabel(ax2)

ax3 = fig.add_subplot(223)
ax3.plot(wsum_in1_to_hid1[:,0], wsum_in1_to_hid1[:,1], label='to hidden 1')
ax3.plot(wsum_in1_to_hid2[:,0], wsum_in1_to_hid2[:,1], label='to hidden 2')
ax3.set_title('Input 1 to hidden layer')
ax3.set_xlabel('Time [s]')
ax3.set_ylabel('Average weight')
#ax3.set_xticks([0]+tics)
ax3.legend()

ax4 = fig.add_subplot(224)
ax4.plot(wsum_in2_to_hid1[:,0], wsum_in2_to_hid1[:,1], label='to hidden 1')
ax4.plot(wsum_in2_to_hid2[:,0], wsum_in2_to_hid2[:,1], label='to hidden 2')
ax4.set_title('Input 2 to hidden layer')
ax4.set_xlabel('Time [s]')
#ax4.set_xticks([0]+tics)
ax4.legend()

plt.tight_layout()
plt.savefig('../../results/xor/Wsum.pdf')
plt.close()



# publication-ready plot
lim_seg = during_learning
r_out = output_rates[lim_seg[0]:lim_seg[1], :]
r_hid1 = hidden1_rates[lim_seg[0]:lim_seg[1], :]
r_hid2 = hidden2_rates[lim_seg[0]:lim_seg[1], :]
r_inp1 = input1_rates[lim_seg[0]:lim_seg[1], :]
r_inp2 = input2_rates[lim_seg[0]:lim_seg[1], :]
er_est = er_trace[lim_seg[0]:lim_seg[1], :]
br_est = br_trace[lim_seg[0]:lim_seg[1], :]

xtik = list(lim_seg[0] + args.durex*np.arange(5)-1)

fig = plt.figure(3, figsize=(7, 7/1.6))
# first subplot is empty

# output layer during learning
ax2 = fig.add_subplot(332)
ax2.plot(r_out[:, 0], r_out[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax2.plot(r_out[:, 0], r_out[:, 1], color=custom_colors['orange'], lw=1, label='BR')
ax2.set_ylabel('Output rates [Hz]')
ax2.set_ylim([0,10])
ax2.set_xticks(xtik)
ax2.spines['top'].set_visible(False)
ax2.legend(loc='best')

ax2_tw = ax2.twinx()
ax2_tw.plot(r_out[:, 0], 100*r_out[:, 1]/r_out[:, 2], color=custom_colors['red'], lw=1)
ax2_tw.plot(er_est[:, 0], 100*br_est[:, 1]/er_est[:, 1], '--', color=custom_colors['red'], lw=1)
ax2_tw.tick_params('y', labelcolor=custom_colors['red'])
ax2_tw.set_yticks([0, 30, 60, 90])
ax2_tw.set_ylabel('Postsyn. BP [\%]', color=custom_colors['red'])
remove_xticklabel(ax2)

#  hidden layer 1
ax3 = fig.add_subplot(334)
ax3.plot(r_hid1[:, 0], r_hid1[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax3.plot(r_hid1[:, 0], r_hid1[:, 1], color=custom_colors['orange'], lw=1, label='BR')
ax3.set_ylabel('Hidden1 rates [Hz]')
ax3.set_xticks((xtik))
ax3.set_ylim([0,10])
ax3.spines['top'].set_visible(False)

ax3_tw = ax3.twinx()
ax3_tw.plot(r_hid1[:, 0], 100*r_hid1[:, 1]/r_hid1[:, 2], color=custom_colors['red'], lw=1)
ax3_tw.tick_params('y', labelcolor=custom_colors['red'])
ax3_tw.set_yticks([0, 30, 60, 90])
ax3_tw.set_ylabel('Postsyn. BP [\%]', color=custom_colors['red'])
remove_xticklabel(ax3)

#  hidden layer 2
ax4 = fig.add_subplot(335)
ax4.plot(r_hid2[:, 0], r_hid2[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax4.plot(r_hid2[:, 0], r_hid2[:, 1], color=custom_colors['orange'], lw=1, label='BR')
ax4.set_ylabel('Hidden2 rates [Hz]')
ax4.set_ylim([0,10])
ax4.set_xticks((xtik))
ax4.spines['top'].set_visible(False)

ax4_tw = ax4.twinx()
ax4_tw.plot(r_hid2[:, 0], 100*r_hid2[:, 1]/r_hid2[:, 2], color=custom_colors['red'], lw=1)
ax4_tw.tick_params('y', labelcolor=custom_colors['red'])
ax4_tw.set_yticks([0, 30, 60, 90])
ax4_tw.set_ylabel('Postsyn. BP [\%]', color=custom_colors['red'])
remove_xticklabel(ax4)


#  input layer 1
ax5 = fig.add_subplot(337)
ax5.plot(r_inp1[:, 0], r_inp1[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax5.set_xlabel("Time [s]")
ax5.set_ylabel('Input1 ER [Hz]')
ax5.set_ylim([0,10])
ax5.set_xticks((xtik))

#  input layer 2
ax6 = fig.add_subplot(338)
ax6.plot(r_inp2[:, 0], r_inp2[:, 2], color=custom_colors['blue'], lw=1, label='ER')
ax6.set_xlabel("Time [s]")
ax6.set_ylabel('Input2 ER [Hz]')
ax6.set_ylim([0,10])
ax6.set_xticks((xtik))


ax7 = fig.add_subplot(333)
lim_seg = before_learning
r_out = output_rates[lim_seg[0]:lim_seg[1]+1, :]
xtik = list(lim_seg[0] + args.durex * np.arange(5)-1)
ax7.plot(r_out[:, 0], r_out[:, 2], color=custom_colors['blue'], lw=1)
ax7.set_ylabel('Output ER [Hz]')
ax7.set_ylim([0,10])
ax7.set_xticks((xtik))


ax8 = fig.add_subplot(336)
lim_seg = after_learning
r_out = output_rates[lim_seg[0]:lim_seg[1]+1, :]
xtik = list(lim_seg[0] + args.durex * np.arange(5)-1)
ax8.plot(r_out[:, 0], r_out[:, 2], color=custom_colors['blue'], lw=1)
ax8.set_xlabel("Time [s]")
ax8.set_ylabel('Output ER [Hz]')
ax8.set_ylim([0,10])
ax8.set_xticks((xtik))


ax9 = fig.add_subplot(339)
ax9.plot(cost_epoch[:,0], cost_epoch[:,1], 'k')
ax9.spines['top'].set_visible(False)
ax9.spines['right'].set_visible(False)
ax9.set_xlabel('Epoch')
ax9.set_ylabel('Cost')


plt.tight_layout()
plt.savefig('../../results/xor/XOR.pdf')
plt.close()


