import numpy as np
import argparse

# parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-numberofneurons", type=int, help="number of neurons", default="1000")
parser.add_argument("-weight", type=float, help="inital weight", default="0.05")
args = parser.parse_args()

N = args.numberofneurons
w = args.weight*np.ones((N,1))

row_indices = np.arange(1,N+1)
row_indices = row_indices[:,np.newaxis]
row_indices = np.broadcast_to(row_indices, (N,2))
X = np.concatenate((row_indices, w), axis=1) 
X = np.concatenate((np.array([[N,N,N]]),X), axis=0)

np.savetxt("../../data/xor/ident_conn.wmat", X, fmt='%d %d %f', header='%%MatrixMarket matrix coordinate real general', comments='')
