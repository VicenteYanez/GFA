#! /usr/bin/python3
"""
@author: Vicente Yáñez
@date: Primer semestre 2016
Universidad de Concepción

Script que corre para crear la lista de estaciones
con sus datos
"""
import sys
import os
import numpy as np
import pdb

from TimeSeriesControl import TimeSeriesControl

codigo = str(sys.argv[1])
lon_min = float(sys.argv[2])
lon_max = float(sys.argv[3])
lat_min = float(sys.argv[4])
lat_max = float(sys.argv[5])
tmin = float(sys.argv[6])
tmax = float(sys.argv[7])

"""
# transformacion de formato intervalo
t1 = fun_tools.string2date(tmin)
t2 = fun_tools.string2date(tmax)
t1 = fun_tools.toYearFraction(t1.year, t1.month, t1.day)
t2 = fun_tools.toYearFraction(t2.year, t2.month, t2.day)
"""

# path archivos
dir_path = os.path.dirname(os.path.realpath(__file__))
lista_gps = '{}/../data/SeriesTiempo/pos_jc'.format(dir_path)
series_dir = '{}/../data/SeriesTiempo/JC/'.format(dir_path)
save_dir = '{}/../example_files/'.format(dir_path)

ts = TimeSeriesControl(codigo, lista_gps, series_dir, save_dir)
ts.load_estations(lon_min, lon_max, lat_min, lat_max, tmin, tmax)
ts.savedata()
