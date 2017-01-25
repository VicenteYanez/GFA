#! /usr/bin/python3

"""
@author: Vicente Yáñez

Generates a time series figure with the gnss and model data
"""
import sys
import os

import numpy as np

from loadGPS import load1stations
from plot_functions import ts_plot

# parameters
codigo = sys.argv[1]
estation = sys.argv[2]

# select directory
home = os.path.dirname(os.path.realpath(__file__))
directory = '{}/../example_files/{}/'.format(home, codigo)

# load and plot functions
tsdata, modeldata = load1stations(estation, directory, True, True)
ts_plot(estation, tsdata[0], tsdata[1:4], modeldata[1:4], savedir=directory)
