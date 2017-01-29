#! /usr/bin/python3

import copy
import os
import pdb

import numpy as np

import fun_vector
from ModeloTrayectoria import ModeloTrayectoria

"""
Functions that operates over a unique model file
"""


def build_model(data, param):
    """
    Genera un nuevo modelo a partir de los parametros dados por el usuario

    Input:
    Lee parametros nuevos desde tmp/cambio_param.txt. Con el mismo
    formato que el guardado por load_stations.
    Codigo fecha

    """
    estacion = data[0]
    serie = data[1]
    # carga de parametros
    modelo = ModeloTrayectoria(serie[0])
    modelo.n = param['polinomio']
    modelo.tjump = np.array(param['saltos'])
    modelo.fperiods = np.array(param['Periodos Fourier'])
    modelo.tlt = np.array(param['Inicio log'])
    modelo.tsc = np.array(param['Escala curva log'])

    # copia del objeto modelo para cada uno de los ejes
    modelo_e = copy.deepcopy(modelo)
    modelo_n = copy.deepcopy(modelo)
    modelo_z = copy.deepcopy(modelo)

    # calculo del modelo de trayectoria
    resultado_e, err_e = modelo_e.modelo_trayectoria(serie[1])
    resultado_n, err_n = modelo_n.modelo_trayectoria(serie[2])
    resultado_z, err_z = modelo_z.modelo_trayectoria(serie[3])

    # resultado total
    res_total_e = (resultado_e[0] + resultado_e[1] +
                   resultado_e[2] + resultado_e[3])
    res_total_n = (resultado_n[0] + resultado_n[1] +
                   resultado_n[2] + resultado_n[3])
    res_total_z = (resultado_z[0] + resultado_z[1] +
                   resultado_z[2] + resultado_z[3])

    # formatear array para guardado
    res_total_e = res_total_e[:, 0]
    res_total_n = res_total_n[:, 0]
    res_total_z = res_total_z[:, 0]

    # ####################################################################
    # GUARDAR MODELO
    # ####################################################################
    resultado_modelo = np.array([serie[0], res_total_e, res_total_n,
                                 res_total_z]).T
    residual = [err_e, err_n, err_z]

    return resultado_modelo, residual


def calc_vector(estacion, file_modelo, vector_file, tipo_vector, *t0):
    """
    Metodo que llama a una funcion para crear un vector
    dependiendo del tipo de calculo elegido
    Input
    *t0: ubicacion de la recta tangente
    """
    modelo = np.loadtxt(file_modelo).T

    # ####################################################################
    # PEDIR VECTOR
    # ####################################################################
    if tipo_vector == 'tangente' and t0:
        t0 = float(t0[0])
        vector, c = fun_vector.tangente(t0, modelo[0], modelo[1],
                                        modelo[2], modelo[3])
    elif tipo_vector == 'ajuste':
        intervalo = [modelo[0][0], modelo[0][-1]]
        vector, c, err = fun_vector.aj_lineal(intervalo, modelo[0],
                                              modelo[1], modelo[2],
                                              modelo[3])
    elif tipo_vector == 'trending':
        print("no implementado")
    else:
        print("tipo de vector no seleccionado")
        return

    # si vector está vacio retornar false
    if vector is False:
        return False

    # GUARDAR VECTOR ######################################################
    # incorporar nombre de estacion a array para guardar
    save_vector = np.hstack(([estacion], vector))

    # guardar vector estacion
    if os.path.isfile(vector_file):
        vectores = np.loadtxt(vector_file, dtype='S5')
        vectores = np.vstack((vectores, save_vector))
        np.savetxt(vector_file, vectores, fmt='%s', delimiter='    ')
    else:
        np.savetxt(vector_file, [save_vector], fmt='%s', delimiter='    ')

    return
