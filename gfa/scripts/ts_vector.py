#! /usr/bin/python3

"""
@author: Vicente Yáñez
@date: Primer semestre 2016
Universidad de Concepción
"""
import sys
import os

import numpy as np

from gfa.gnss_analysis.model_accions import calc_vector
from gfa.load_param import Config


def main(identifier, estacion, vector_type, t):
    # path archivos
    save_dir = '{}{}/'.format(Config.config['PATH']['output_dir'],
                              identifier)
    vector_file = '{}vectors.txt'.format(save_dir)
    file_modelo = '{}modelo/{}.txt'.format(save_dir, estacion)

    if vector_type == "tangent":
        t0 = t[0]
        res = calc_vector(estacion, file_modelo, vector_file, vector_type, t0)
    elif vector_type == "fit":
        res = calc_vector(estacion, file_modelo, vector_file, vector_type, t)
