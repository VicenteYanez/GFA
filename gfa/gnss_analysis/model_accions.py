#! /usr/bin/python3

import copy
import os
import json

import numpy as np
import pandas as pd

from gfa.gnss_analysis import fun_vector
from gfa.gnss_analysis.ModeloTrayectoria import ModeloTrayectoria

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
    res_total_e = modelo_e.total
    res_total_n = modelo_n.total
    res_total_z = modelo_z.total

    # ####################################################################
    # GUARDAR MODELO
    # ####################################################################
    modelresult = np.array([serie[0], res_total_e, res_total_n, res_total_z,
                            resultado_e[0], resultado_n[0], resultado_z[0],
                            resultado_e[1], resultado_n[1], resultado_z[1],
                            resultado_e[2], resultado_n[2], resultado_z[2],
                            resultado_e[3], resultado_n[3], resultado_z[3]]).T
    residual = [err_e, err_n, err_z]

    return modelresult, residual


def calc_vector(estacion, file_modelo, vector_file, vector_type, aux=False):
    """
    Method that calls call the functions for compute
    a velocity vector, depending of the . And finally saves the
    the vector in vectors.txt file.
    Input
    aux: It is a value that depends of the vector_type variable.
        If it is "fit":
            Time range for the lineal fit (list of lenght 2)
        If it is "tangent":
            location of the tangent line (float)
    """
    modelo = np.loadtxt(file_modelo).T

    # ####################################################################
    # PEDIR VECTOR
    # ####################################################################
    if vector_type == 'tangent':
        t0 = float(aux)

        # Si la serie no contiene el tiempo t, regresa falso
        if modelo[0][len(modelo[0])-1] < float(t0):
            print('Error, station do not have the time {}'.format(t0))
            return False

        vector, c = fun_vector.tangente(t0, modelo[0], modelo[1],
                                        modelo[2], modelo[3])
        t1 = t0
        t2 = t0
    elif vector_type == 'fit':
        trange = aux

        # Si la serie no contiene el tiempo t, regresa falso
        if modelo[0][len(modelo[0])-1] < float(trange[0]):
            print('Error, station do not have the time {}'.format(trange[0]))
            return False

        vector, c, err = fun_vector.fit(trange, modelo[0], modelo[1],
                                        modelo[2], modelo[3])
        t1 = trange[0]
        t2 = trange[1]
    elif vector_type == 'trending':
        print("not implemented yet")
    else:
        print("vector type don't selected")
        return False

    # si vector está vacio retornar false
    if vector is False:
        return False

    # SAVE VECTOR ######################################################
    # create dataframe with the new data
    savedf = pd.DataFrame({'station': [estacion], 'vector_type': [vector_type],
                           'vector_e': [vector[0]], 'vector_n': [vector[1]],
                           'vector_z': [vector[2]], 'c_e': [c[0]],
                           'c_n': [c[1]], 'c_z': [c[2]],
                           'start_time': [t1], 'end_time': [t2]})
    # order columns
    savedf = savedf[['station', 'vector_type', 'vector_e', 'vector_n',
                     'vector_z', 'c_e', 'c_n', 'c_z',
                     'start_time', 'end_time']]
    # add it with the previus data
    if os.path.isfile(vector_file):
        old_data = pd.read_csv(vector_file)
        old_data = old_data[['station', 'vector_type', 'vector_e', 'vector_n',
                             'vector_z', 'c_e', 'c_n', 'c_z',
                             'start_time', 'end_time']]
        new_data = old_data.append(savedf)
        new_data.to_csv(vector_file, index=False)

    else:
        # if there is no data saved previusly
        savedf.to_csv(vector_file, index=False)

    return True


def save_model(modelo, model_file):
    """
    Save the data of one model
    """
    head = 'Time[yr], Est[mm], North[mm], Vertical[mm],\
Poly_e, Poly_n, Poly_z, Jump_e, Jump_n, Jump_z, Fourier_e, Fourier_n,\
Fourier_z, Log_e, Log_n, Log_z'
    # save = np.array(modelo[1]).T
    np.savetxt(model_file, modelo, fmt='%s', header=head, delimiter='    ')
    return


def upgrade_list(estacion, parametros, residual, directory):
    """
    Upgrade the model list with the new parameters
    """
    # pasar residual y parametros to json
    residual = {'Residual E': residual[0].tolist(),
                'Residual N': residual[1].tolist(),
                'Residual Z': residual[2].tolist()}
    parametros_json = json.dumps(parametros)
    residual_json = json.dumps(residual)
    head = 'station    longitude    latitude    Parameters    Error'
    archivo = '{}modelo_lista.txt'.format(directory)
    data = np.loadtxt(archivo, dtype=bytes, delimiter='    ').astype(str)
    for i, row in enumerate(data):
        if row[0] == estacion:
            data[i] = np.array([estacion, row[1], row[2],
                                parametros_json, residual_json])
    # actualizar txt
    np.savetxt(archivo, data, fmt='%s', delimiter='    ', header=head)
    return


def load_vector(vectorfile, station):
    """
    Function that load the vector from the vector.txt file
    """
    estlist = np.loadtxt(vectorfile, usecols=[0], delimiter=',', dtype=bytes,
                         skiprows=1, ndmin=1).astype(str)
    vectors = np.loadtxt(vectorfile, usecols=[2, 3, 4], delimiter=',',
                         skiprows=1, ndmin=2, dtype=float)
    c = np.loadtxt(vectorfile, usecols=[5, 6, 7], skiprows=1, dtype=float,
                   ndmin=2, delimiter=',')
    t1t2 = np.loadtxt(vectorfile, usecols=[8, 9], skiprows=1, dtype=float,
                      ndmin=2, delimiter=',')
    data = []
    try:
        for i, est in enumerate(estlist):
            if est == station:
                # list with t inicial and t final
                # if vector=tangent, t inicial == t final
                tmin_tmax = [t1t2[i][0], t1t2[i][-1]]
                data.append([vectors[i], c[i], tmin_tmax])

    except TypeError:
        print('there is no data in {}'.format(vectorfile))
        exit()
    return data
