#! /usr/bin/python3
"""
@author: Vicente Yáñez
@date: 2017

Extract the time series data giving a certain time range and area

"""
import sys
import os
import numpy as np
import pdb
import shutil
import logging

from TimeSeriesControl import TimeSeriesControl
from ModelControl import ModelControl

# login config
logging.basicConfig(filename='../gfa.log', level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(name)s %(message)s')
logger = logging.getLogger(__name__)

try:
    codigo = str(sys.argv[1])
    lon_min = float(sys.argv[2])
    lon_max = float(sys.argv[3])
    lat_min = float(sys.argv[4])
    lat_max = float(sys.argv[5])
    tmin = float(sys.argv[6])
    tmax = float(sys.argv[7])

    if lon_min > lon_max or lat_min > lat_max or tmin > tmax:
        raise ValueError()

except(ValueError, IndexError) as err:
    print('A wild error had raised when GFA was reading the parameters!\
 Please, check your input parameter and restart the script')
    logger.error(err)
    exit()

# path archivos
dir_path = os.path.dirname(os.path.realpath(__file__))
lista_gps = '{}/../data/TimeSeries/station_list.txt'.format(dir_path)
series_dir = '{}/../data/TimeSeries/txtfiles/'.format(dir_path)
save_dir = '{}/../example_files/'.format(dir_path)
# optional: dir of the trajectory model
model_dir = '{}/../example_files/general_solution/Modelo/'.format(dir_path)
m_ls = '{}/../example_files/general_solution/resume.txt'.format(dir_path)

# check if save_dir exist
select_dir = '{}{}'.format(save_dir, codigo)
if os.path.exists(select_dir):
    erase = input('{} exist. Do you want to erase it? y/n(n)'.format(
        select_dir))
    if erase == 'y':
        shutil.rmtree(select_dir)
        print('Directory removed')
    else:
        rename = input('Rename the identificator parameter:')
        codigo = rename
try:
    ts = TimeSeriesControl(codigo, lista_gps, series_dir, save_dir)
    ts.load_estations(lon_min, lon_max, lat_min, lat_max, tmin, tmax)
    ts.savedata()

    # optional: if you have a trajectory model you can load the same way
    if os.path.exists(model_dir):
        model = ModelControl(codigo, m_ls, model_dir, save_dir)
        model.load_model(lon_min, lon_max, lat_min, lat_max, tmin, tmax)
        model.savedata()

except:
    print('A wild error had raised! Please, check your input param or go and\
see the log')
    logger.error(err)
