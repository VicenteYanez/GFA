#! /usr/bin/python3

"""
@author: Vicente Yáñez

Generates a time series figure with the gnss and model data
Parameters
[1] identifier code
[2] estation to plot
[3] If the third parameter is -m, besides the time series, the model
is plot.
[4] If the forth parameter is -v, the velocity vector is plot
"""

import sys
import os
import pdb

import matplotlib.pyplot as plt
import numpy as np

from loadGPS import load1stations
from model_accions import load_vector
import fun_tsplot as pltfun
import log_config as log


# parameters
codigo = sys.argv[1]
estation = sys.argv[2]

# select directory
home = os.path.dirname(os.path.realpath(__file__))
directory = '{}/../example_files/{}/'.format(home, codigo)
savefile = '{}{}.png'.format(directory, estation)

# load data and plot function
tsdata, modeldata = load1stations(estation, directory, True, True)
f, axes = pltfun.plot(estation, tsdata[1][0], tsdata[1][1:4])

# optional plots: add model
try:
    if sys.argv[3] == '-m':
        f, axes = pltfun.add_modelo(f, axes, modeldata[1][0],
                                    modeldata[1][1:4])
        print('Model added')
except IndexError:
    print('Model not added')
    pass

# optional plots: add vector
try:
    if sys.argv[4] == '-v':
        # file vector
        filevector = '{}vectors.txt'.format(directory)
        try:
            # load vectors
            vectors = load_vector(filevector, estation)
            # plot
            f, axes = pltfun.add_velocity(f, axes, tsdata[1][0], vectors)
            print('Vector added')
        except FileNotFoundError as err:
            print("Error: Vector file not found")
            log.logger.error(err)

except IndexError:
    print('Vector not added')
    pass

# save plot
pltfun.save_figure(f, savefile)
