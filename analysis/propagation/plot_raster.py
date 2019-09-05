import matplotlib.pyplot as plt
import sys
sys.path.append('../')
import raster_analysis as ra
import voltage_analysis as va
import numpy as np
import argparse
import utils_plotting as up

parser = argparse.ArgumentParser()
parser.add_argument("-datadir", type=str, help="data directory", default="")
parser.add_argument("-resultdir", type=str, help="parent result directory", default="")

args = parser.parse_args()

rasterfile_pop1 = args.datadir + 'noise_matching.0.ras1_seed1'
rasterfile_pop2 = args.datadir + 'noise_matching.0.ras2_seed1'

ra.display_raster_with_eventplot(rasterfile_pop1, args.resultdir + 'Raster1.pdf')
ra.display_raster_with_eventplot(rasterfile_pop2, args.resultdir + 'Raster2.pdf')


