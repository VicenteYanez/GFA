#! /usr/bin/python3
"""
@author: Vicente Yáñez
@date: 2016
Universidad de Concepción

Calculate the trajectory model for all the stations using the same
parameters
"""

import os
import numpy as np

from gfa.gnss_analysis.ModelControl import ModelControl
from gfa.load_param import Config
import gfa.log_config as log


def main():
    # Path archivos
    lista_gps = Config.config['PATH']['ListaGPS']
    series_dir = Config.config['PATH']['GPSdata']
    save_dir = Config.config['PATH']['output_dir']
    eq_file = Config.config['PATH']['eqfile']

    # area y tiempo
    lon = [-180, 180]
    lat = [-90, 90]
    t = [1950, 2050]

    model = ModelControl('general_solution', lista_gps, series_dir, save_dir)
    data, lista = model.load_estations(lon[0], lon[1], lat[0], lat[1], t[0],
                                       t[1])
    model.build_model_all(1, [1], eq_file)
