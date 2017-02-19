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

import fun_datos
from loadGPS import load1stations
from model_accions import build_model, save_model, upgrade_list

# parametros
codigo = sys.argv[1]
estation = sys.argv[2]

path = os.path.dirname(os.path.abspath(__file__))
directory = '{}/../example_files/{}/'.format(path, codigo)
model_file = '{}modelo/{}.txt'.format(directory, estation)

if os.path.isfile(model_file) is False:
    sys.exit('Error, model had not been loaded')

# parametros del nuevo modelo
parametros = {"polinomio": 1,
              "saltos": [2010.01, 2014.24832],
              "Escala curva log": [1],
              "Periodos Fourier": [1],
              "Inicio log": [2014.24832]}

tsdata, modeldata = load1stations(estation, directory, True, model=False)
modelo, residual = build_model(tsdata, parametros)
save_model(modelo, model_file)
upgrade_list(estation, parametros, residual, directory)
