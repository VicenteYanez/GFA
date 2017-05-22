#! /usr/bin/python3

"""
@author: Vicente Yáñez

Generates a time series figure with the gnss and model data
Parameters
[1] identifier code
[2] station to plot
[3] If the third parameter is -m, besides the time series, the model
is plot.
[4] If the forth parameter is -v, the velocity vector is plot
"""

import sys
import os
import pdb

import matplotlib.pyplot as plt
import numpy as np

from gfa.gnss_analysis.loadGPS import load1stations
from gfa.gnss_analysis.model_accions import load_vector
import gfa.gnss_analysis.fun_tsplot as pltfun
import gfa.log_config as log
from gfa.load_param import Config


def main(codigo, station, addmodel, addvector):
    # select directory
    directory = '{}{}/'.format(Config.config['PATH']['output_dir'], codigo)
    savefile = '{}{}.png'.format(directory, station)

    # load data and plot function
    tsdata, modeldata = load1stations(station, directory, True, True)
    f, axes = pltfun.plot(station, tsdata[1][0], tsdata[1][1:4])

    # optional plots: add model
    if addmodel == '-m':
        f, axes = pltfun.add_modelo(f, axes, modeldata[1][0],
                                    modeldata[1][1:4])
        print('Model added')
    else:
        print('Model not added')

    # optional plots: add vector
    if addvector == '-v':
        # file vector
        filevector = '{}vectors.txt'.format(directory)
        try:
            # load vectors
            vectors = load_vector(filevector, station)
            # plot
            f, axes = pltfun.add_velocity(f, axes, tsdata[1][0], vectors)
            print('Vector added')
        except FileNotFoundError as err:
            print("Error: Vector file not found")
            log.logger.error(err)

    else:
        print('Vector not added')

    # save plot
    pltfun.save_figure(f, savefile)
