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

from ModelControl import ModelControl


# Path archivos
dir_path = os.path.dirname(os.path.realpath(__file__))
lista_gps = '{}/../data/SeriesTiempo/pos_jc'.format(dir_path)
series_dir = '{}/../data/SeriesTiempo/JC/'.format(dir_path)
save_dir = '{}/../example_files/'.format(dir_path)
eq_file = '{}/../data/eq_file.txt'.format(dir_path)

# area y tiempo
lon = [-180, 180]
lat = [-90, 90]
t = [1950, 2050]

model = ModelControl('general_solution', lista_gps, series_dir, save_dir)
data, lista = model.load_estations(lon[0], lon[1], lat[0], lat[1], t[0], t[1])
model.build_model_all(1, [1], eq_file)
