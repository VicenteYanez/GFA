#! /usr/bin/python3

"""
@author: Vicente Yáñez

Script que se ejecuta cada vez que se necesite modificar el modelo
de una estación
"""

import sys
import os
import pdb

import numpy as np

import gfa.gnss_analysis.fun_datos
from gfa.gnss_analysis.loadGPS import load1stations
from gfa.load_param import Config
import gfa.gnss_analysis.model_accions as model_accions


def main(codigo, station, polinomio, steps, fourier, logscale, logstart):
    directory = '{}{}/'.format(Config.config['PATH']['output_dir'], codigo)
    model_file = '{}modelo/{}.txt'.format(directory, station)

    if os.path.isfile(model_file) is False:
        sys.exit('Error, model had not been loaded')

    # parametros del nuevo modelo
    parametros = {"polinomio": polinomio,
                  "saltos": steps,
                  "Escala curva log": logscale,
                  "Periodos Fourier": fourier,
                  "Inicio log": logstart}

    tsdata, modeldata = load1stations(station, directory, True, model=False)
    modelo, residual = model_accions.build_model(tsdata, parametros)
    model_accions.save_model(modelo, model_file)
    model_accions.upgrade_list(station, parametros, residual, directory)
