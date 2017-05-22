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


from gfa.gnss_analysis.TimeSeriesControl import TimeSeriesControl
from gfa.gnss_analysis.ModelControl import ModelControl
from gfa.load_param import Config
import gfa.log_config as log


def main(codigo, lon_min, lon_max, lat_min, lat_max, tmin, tmax):
    try:
        if lon_min > lon_max or lat_min > lat_max or tmin > tmax:
            raise ValueError()

    except(ValueError, IndexError) as err:
        print('A wild error had raised when GFA was reading the parameters!\
     Please, check your input parameter and restart the script')
        log.logger.error(err)
        exit()

    # path archivos
    lista_gps = Config.config['PATH']['ListaGPS']
    series_dir = Config.config['PATH']['GPSdata']
    save_dir = Config.config['PATH']['output_dir']
    # optional: dir of the trajectory model
    model_dir = '{}general_solution/Modelo/'.format(save_dir)
    m_ls = '{}general_solution/resume.txt'.format(save_dir)

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

    except() as err:
        print('A wild error had raised! Please, check your input param or go\
              and see the log')
        log.logger.error(err)
