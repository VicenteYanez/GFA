#! /usr/bin/python3

"""
@author: Vicente Yáñez
@date: Diciembre 2016

Set of functions to load and modify GNSS time series
"""

import os
import numpy as np
import pdb


def clean_serie(serie, err_max=10):
    """
    Funcion que limpia las mediciones dentro de una serie de tiempo
    con un error mayor al predefinido
    """
    serie_limpia = []
    for i, med in enumerate(serie):
        # transformacion de string a float
        err_e = serie[i][4]
        err_n = serie[i][5]
        err_z = serie[i][6]
        if (err_e < err_max and err_n < err_max and err_z < err_max):
            serie_limpia.append(serie[i])
    return serie_limpia


def select_time(serie, intervalo):
    """
    Función que corta la serie para un intervalo definido
    Input:
    serie     : serie de tiempo array-like format_csn
    intervalo : intervalo de tiempo pedido array[min, max]
    """
    serie2 = []
    for i, med in enumerate(serie):
        t_serie = serie[i][0]
        if t_serie > intervalo[0] and t_serie < intervalo[1]:
            serie2.append(serie[i])
    return np.array(serie2)


def format_csn(lon, lat, intervalo, estac_list, path_series):
    """
    Function that load all the time series in a directory. Returns one list
    with the data and other list with the name of stations and his locations.
    The method can only load the data if they are in the format giving for
    the Centro Sismologico Nacional (CSN).
    ### Input
    * lon: [list] with the min and maximum value of longitude.
    * lat: [list] with the min and maximun value of latitude.
    * intervalo: with the min and maximun value of time (in years).
    * estac_list: [str] with the path to a text file with the name and
    geographic coordinates of each stations.
    * path_series [str] with the path to a directory that have the
    time series's file.
    ### Output
    * data = lenght x nº of stations

                    data[0] : ['name estation', [0, ..., 6]]
                    Where 0 to 6:
                        0. time
                        1. East desplacement
                        2. North desplacement
                        3. Vertical desplacement
                        4. East Error
                        5. North Error
                        6. Vertical Error
                    example:
                        list[4][1][2] : north desplacment of station number 4
                        list[4][0]    : name of station number 4
    * list_estations: [list] with the name, longitude and latitude
    for each station.
    """
    # Carga de lista de estaciones
    estaciones = np.loadtxt(estac_list, usecols=[0], dtype='S5', skiprows=1)
    lon_est = np.loadtxt(estac_list, usecols=[1], skiprows=1, dtype=float)
    lat_est = np.loadtxt(estac_list, usecols=[2], skiprows=1, dtype=float)

    # array para guardar datos
    data = np.array([])
    lista_estac = []
    data_vacio = True

    # ###################################################################
    # CICLO DE CARGA
    # ###################################################################
    for i, estacion in enumerate(estaciones):
        estacion = estacion.decode('utf-8')
        # ###################################################
        # Seleccion por area
        if lat_est[i] < lat[0] or lat_est[i] > lat[1]:
            continue
        if lon_est[i] < lon[0] or lon_est[i] > lon[1]:
            continue
        # ###################################################

        # path estacion i, para cada eje
        dir_estac = "{}{}.txt".format(path_series, estacion)

        # verifica que exista la estacion
        if ((os.path.isfile(dir_estac) is False)):
            print("Faltan archivos de {}".format(estacion))
            continue
        # verifica que estacion i tenga un tamaño minimo
        stinfo = os.stat(dir_estac)
        if stinfo.st_size > 60:
            # abrir archivos de estacion
            data_estac = np.loadtxt(dir_estac,
                                    usecols=(1, 2, 3, 4, 5, 6, 7, 8),
                                    dtype=float)
            # carga de datos
            year = data_estac.T[0]
            dias = data_estac.T[1] / 365.25
            tiempo = year + dias
            desp_e = data_estac.T[2]
            desp_n = data_estac.T[3]
            desp_u = data_estac.T[4]
            desp_e_err = data_estac.T[5]
            desp_n_err = data_estac.T[6]
            desp_u_err = data_estac.T[7]

        # nota: python pasa los valores a string en esta operacion
        # corregir usando .astype(np.float).T
        data_temp = np.array([tiempo, desp_e, desp_n, desp_u, desp_e_err,
                              desp_n_err, desp_u_err]).T
        # limpiar datos
        data_temp_limpia = clean_serie(data_temp)
        data_temp2 = select_time(data_temp_limpia, intervalo)

        # si se tienen menos de 10 mediciones descarta la estacion
        if len(data_temp2) < 10:
            print("Pocos datos para estacion {}".format(estacion))
            continue

        data_temp2 = list(data_temp2.T)

        lista_estac.append([estacion, lon_est[i], lat_est[i]])

        # la primera vez data = data_temp, luego va apilando
        if data_vacio:
            data = [[estacion, data_temp2]]
            data_vacio = False
        else:
            data.append([estacion, data_temp2])

    return data, lista_estac


def format_model(lon, lat, intervalo, estac_list, path_series):
    """
    Similar to format_csn but load the trajectory model instead.
    ### Input
    * lon: [list] with the min and maximum value of longitude.
    * lat: [list] with the min and maximun value of latitude.
    * intervalo: with the min and maximun value of time (in years).
    * estac_list: [str] with the path to a text file with the name and
    geographic coordinates of each stations.
    * path_series [str] with the path to a directory that have the
    time series's file.
    ### Output
    * data = lenght x nº of stations

                    data[0] : ['name estation', [0, ..., 6]]
                    Where 0 to 6:
                        0. time
                        1. East desplacement
                        2. North desplacement
                        3. Vertical desplacement
                    example:
                        data[4][1][2] : north desplacment of station number 4
                        data[4][0]    : name of station number 4
    * list_estations: [list] with the name, longitude and latitude,
                      parameters and the residual of the model
                      for each station.
    """
    # Carga de lista de estaciones
    estaciones = np.loadtxt(estac_list, usecols=[0], dtype=bytes,
                            skiprows=1).astype(str)
    lon_est = np.loadtxt(estac_list, usecols=[1], skiprows=1, dtype=float)
    lat_est = np.loadtxt(estac_list, usecols=[2], skiprows=1, dtype=float)
    param = np.loadtxt(estac_list, usecols=[3], skiprows=1, dtype=bytes,
                       delimiter='    ').astype(str)
    residual = np.loadtxt(estac_list, usecols=[4], skiprows=1, dtype=bytes,
                          delimiter='    ').astype(str)

    # array para guardar datos
    data = np.array([])
    lista_estac = []
    data_vacio = True

    # ###################################################################
    # CICLO DE CARGA
    # ###################################################################
    for i, estacion in enumerate(estaciones):
        # ###################################################
        # Seleccion por area
        if lat_est[i] < lat[0] or lat_est[i] > lat[1]:
            continue
        if lon_est[i] < lon[0] or lon_est[i] > lon[1]:
            continue
        # ###################################################

        # path estacion i, para cada eje
        dir_estac = "{}{}.txt".format(path_series, estacion)

        # verifica que exista la estacion
        if ((os.path.isfile(dir_estac) is False)):
            print("Faltan archivos de {}".format(estacion))
            continue
        # verifica que estacion i tenga un tamaño minimo
        stinfo = os.stat(dir_estac)
        if stinfo.st_size > 60:
            # abrir archivos de estacion
            data_estac = np.loadtxt(dir_estac, usecols=(0, 1, 2, 3),
                                    dtype=float, skiprows=1)
            # carga de datos
            tiempo = data_estac.T[0]
            desp_e = data_estac.T[1]
            desp_n = data_estac.T[2]
            desp_u = data_estac.T[3]

        data_temp = np.array([tiempo, desp_e, desp_n, desp_u]).T
        # select by time
        data_temp2 = select_time(data_temp, intervalo)

        # si se tienen menos de 10 mediciones descarta la estacion
        if len(data_temp2) < 10:
            print("Pocos datos para estacion {}".format(estacion))
            continue

        data_temp2 = list(data_temp2.T)
        lista_estac.append([estacion, lon_est[i], lat_est[i], param[i],
                            residual[i]])

        # la primera vez data = data_temp, luego va apilando
        if data_vacio:
            data = [[estacion, data_temp2]]
            data_vacio = False
        else:
            data.append([estacion, data_temp2])

    return data, lista_estac


def load1stations(name, directory, ts=True, model=False):
    """
    Function that load the data generated by the TimeSeriesControl class
    for the station 'name'.

    ts    : by defect True. Indicates if the function load the data of the
            time series.
    model : by defect False. Indicates if the function load the data of the
            model.
    """
    if ts is False and model is False:
        print('No data selected, time series and model are false.')
        return None
    tsdata = False
    modeldata = False
    if ts is True:
        ts_dir = '{}series/{}.txt'.format(directory, name)
        tsdata = np.loadtxt(ts_dir, dtype=float, skiprows=1)
        tsdata = [name, tsdata.T]
    if model is True:
        model_dir = '{}modelo/{}.txt'.format(directory, name)
        modeldata = np.loadtxt(model_dir, dtype=float, skiprows=1)
        modeldata = [name, modeldata.T]
    return tsdata, modeldata
