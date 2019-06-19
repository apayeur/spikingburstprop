import sys
sys.path.append('../')
import raster_analysis as ra
import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-filesuffix", type=str, help="suffix indicating feedback weights", default="")
parser.add_argument("-resultdir", type=str, help="parent result directory", default=".")
parser.add_argument("-datadir", type=str, help="parent data directory", default=".")

args = parser.parse_args()


currents = np.arange(-400, 1600., 200.)

# BP, ER, BR for population 1 and 2
filenames1 = [args.datadir + '/propagation_ficurve.0.brate1_somacurrent0.000000',
              args.datadir + '/propagation_ficurve.0.brate1_somacurrent100.000000']
filenames2 = [args.datadir + '/propagation_ficurve.0.brate2_somacurrent0.000000',
              args.datadir + '/propagation_ficurve.0.brate2_somacurrent100.000000']

outfile = args.resultdir + '/BP1versusIdend' + args.filesuffix + '.pdf'
ra.display_responsecurves(currents, outfile, *filenames1)
outfile = args.resultdir + '/BP2versusIdend' + args.filesuffix + '.pdf'
ra.display_responsecurves(currents, outfile, *filenames2)

# PV and SOM rates
'''
filenames_pv = [args.datadir + '/propagation_ficurve.0.pvrate_somacurrent100.000000']
filenames_som = [args.datadir + '/propagation_ficurve.0.somrate_somacurrent100.000000']
ra.display_ficurves(currents, args.resultdir + '/SOMFRversusIdend' + args.filesuffix + '.pdf', *filenames_som)
ra.display_ficurves(currents, args.resultdir + '/PVFRversusIdend' + args.filesuffix + '.pdf', *filenames_pv)
'''