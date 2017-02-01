#! /usr/bin/python3

"""
@author: Vicente Yáñez
@date: Primer semestre 2016
Universidad de Concepción

Script que se ejecuta para calcular el vector de
desplazamiento para una estacion en un periodo de tiempo

Argumentos
[1] : codigo de almacenado de archivo
[2] : nombre estacion
[3] : tipo de vector elegido (tangente, ajuste, trending)
[4] : opcional, float t para vector tangente
"""
import sys
import os

import numpy as np

from model_accions import calc_vector


# Argumentos
identifier = sys.argv[1]
estacion = sys.argv[2]
tipo_vector = sys.argv[3]

# path archivos
dir_path = os.path.dirname(os.path.realpath(__file__))
save_dir = '{}/../example_files/{}/'.format(dir_path, identifier)
vector_file = '{}vectors.txt'.format(save_dir)
file_modelo = '{}modelo/{}.txt'.format(save_dir, estacion)

if len(sys.argv) > 4:
    arg = sys.argv[4]
    res = calc_vector(estacion, file_modelo, vector_file, tipo_vector, arg)
else:
    res = calc_vector(estacion, file_modelo, vector_file, tipo_vector)

if res is False:
    print('ts_vector was failed')
