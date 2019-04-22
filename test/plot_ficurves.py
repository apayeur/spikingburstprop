import sys
sys.path.append('../analysis/')
import raster_analysis as ra
import numpy as np

# brate files
filenames = ['./test-output/ficurve_burstpoisson.0.brate_somacurrent0.000000',
             './test-output/ficurve_burstpoisson.0.brate_somacurrent100.000000',
             './test-output/ficurve_burstpoisson.0.brate_somacurrent200.000000']

currents = np.arange(-200, 1500., 100.)
outfile = './test-output/BPversusIdend.pdf'
ra.display_ficurves(currents, outfile, *filenames)
