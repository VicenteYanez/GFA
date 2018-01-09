#! /usr/bin/python3

"""
@author: Vicente Yáñez

Script que se ejecuta cada vez que se necesite modificar el modelo
de una estación
"""
from datetime import datetime
import time
import sys
import os

from gfa.data_tools.loadGPS import load1stations
from gfa.load_param import Config
import gfa.gnss_analysis.model_accions as model_accions


def toYearFraction(date):
    def sinceEpoch(date):  # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction


def main(codigo, station, polinomio, jumps, fourier, logscale, logstart):
    directory = '{}{}/'.format(Config.config['PATH']['output_dir'], codigo)
    model_file = '{}modelo/{}.txt'.format(directory, station)

    if os.path.isfile(model_file) is False:
        sys.exit('Error, model had not been loaded')

    if fourier == '':
        fourier = []
    else:
        fourier = [float(s) for s in fourier.split(',')]
    if jumps == '':
        jumps1 = []
        jumps2 = []
    else:
        jumps1 = [float(toYearFraction(datetime.strptime(s, "%Y-%m-%d")))
                  for s in jumps.split(',')]
        jumps2 = [s for s in jumps.split(',')]
    logscale = [float(s) for s in logscale.split(',')]
    if logstart == '':
        logstart1 = []
        logstart2 = []
    else:
        logstart1 = [float(toYearFraction(datetime.strptime(s, "%Y-%m-%d")))
                     for s in logstart.split(',')]
        logstart2 = [s for s in logstart.split(',')]

    # parametros del nuevo modelo
    # number 1 is for the trajectory model calculation
    # number 2 is for save the parameters in the .txt file
    parametros1 = {"polinomio": polinomio,
                   "saltos": jumps1,
                   "Escala curva log": logscale,
                   "Periodos Fourier": fourier,
                   "Inicio log": logstart1}

    parametros2 = {"polinomio": polinomio,
                   "saltos": jumps2,
                   "Escala curva log": logscale,
                   "Periodos Fourier": fourier,
                   "Inicio log": logstart2}

    tsdata, modeldata = load1stations(station, directory, True, model=False)
    modelo, residual = model_accions.build_model(tsdata, parametros1)
    model_accions.save_model(modelo, model_file)
    model_accions.upgrade_list(station, parametros2, residual, directory)
