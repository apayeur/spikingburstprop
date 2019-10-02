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
data_file_prefix = '../../data/xor/xor.0.wsum_'
wsum_noise_to_out  = np.loadtxt(data_file_prefix + 'noise_to_out')
wsum_noise_to_hid1 = np.loadtxt(data_file_prefix + 'noise_to_hid1')
wsum_noise_to_hid2 = np.loadtxt(data_file_prefix + 'noise_to_hid2')

after_learning = [wsum_noise_to_out.shape[0]-int(4*args.durex/1.)+1, wsum_noise_to_out.shape[0]+1]
before_learning = [int(3*args.alpha/1.)+1, int(3*args.alpha/1.)+int(4*args.durex/1.)+1]
during_learning = [int(3*args.alpha/1.)+int(args.numex//2*args.durex/1.)+1,
                   int(3*args.alpha/1.)+int(args.numex//2*args.durex/1.)+int(4*args.durex/1.)+1]


fig = plt.figure(figsize=(3, 6./1.6))

plt.subplot(311)
plt.plot(wsum_noise_to_out[:,0], wsum_noise_to_out[:,1], 'k')
plt.ylabel("Bias output layer")
plt.xlabel("Time [s]")

plt.subplot(312)
plt.plot(wsum_noise_to_hid1[:,0], wsum_noise_to_hid1[:,1], 'k')
plt.ylabel("Bias hidden 1")
plt.xlabel("Time [s]")

plt.subplot(313)
plt.plot(wsum_noise_to_hid2[:,0], wsum_noise_to_hid2[:,1], 'k')
plt.ylabel("Bias hidden 2")
plt.xlabel("Time [s]")

plt.tight_layout()
plt.savefig('../../results/xor/Bias.pdf')
plt.close()