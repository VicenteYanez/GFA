#! /usr/bin/python3

"""
@author: Vicente Yáñez
@date: Diciembre 2016

Funciones relacionadas a cargar y modificar series de tiempo GPS
"""

import os
import numpy as np
import pdb


def _limpiar_serie(serie, err_max=10):
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


def _extraer_tiempo(serie, intervalo):
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
    Metodo que carga todas las series de tiempo con el formato del CSN.
    Entrega array con lista de nombres de estaciones.

    format      :   nombre  año dia eje(lon,lat,z), errores (lon, lat, z)

    Input
        lon         :   array-lide [lon min, lon max]
        lat         :   array-like [lat min, lat max]
        intervalo   :   array-like [tiempo min, tiempo max]
        estac_list  :   archivo txt con el nombre y ubicacion de cada estacion
        path_series :   directorio donde estan las series de tiempo

    Output:
    data        :   array de 2d, de longitud num medidas x 7
                    Columna
                        1. Nombre estacion
                        2. tiempo
                        3. Desplazamiento Este
                        4. Desplazamiento Norte
                        5. Desplazamiento Vertical
                        6. Error Este
                        7. Error Norte
                        8. Error Vertical
    """
    # Carga de lista de estaciones
    estaciones = np.loadtxt(estac_list, usecols=[0], dtype='S5', skiprows=1)
    lat_est = np.loadtxt(estac_list, usecols=[2], skiprows=1, dtype=float)
    lon_est = np.loadtxt(estac_list, usecols=[1], skiprows=1, dtype=float)

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
        data_temp_limpia = _limpiar_serie(data_temp)
        data_temp2 = _extraer_tiempo(data_temp_limpia, intervalo)

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


def carga_modelo(lon, lat, intervalo, estac_list, path_modelo):
    """
    Metodo que carga la solucion al modelo de trayectoria
    Entrega array con lista de nombres de estaciones.

    format      :   tiempo  desplazamiento

    Input
        lon         :   array-lide [lon min, lon max]
        lat         :   array-like [lat min, lat max]
        intervalo   :   array-like [tiempo min, tiempo max]
        estac_list  :   archivo txt con el nombre y ubicacion de cada estacion
        path_modelo :   directorio con los archivos del modelo de trayectoria

    Output:
    data        :   array de 2d, de longitud num medidas x 7
                    Columna
                        1. Nombre estacion
                        2. tiempo
                        3. Desplazamiento Este
                        4. Desplazamiento Norte
                        5. Desplazamiento Vertical
    """
    # Carga de lista de estaciones
    estaciones = np.loadtxt(estac_list, usecols=[0],
                            dtype=str, skiprows=1)
    lat_est = np.loadtxt(estac_list, usecols=[2], skiprows=1)
    lon_est = np.loadtxt(estac_list, usecols=[1], skiprows=1)

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
        if lat_est[i] < lat[0] and lat_est[i] > lat[1]:
            continue
        if lon_est[i] < lon[0] and lon_est[i] > lon[1]:
            continue
        # ###################################################

        # path estacion i, para cada eje
        estac_file = path_modelo + estacion
        estac_file_e = estac_file + '_e.txt'
        estac_file_n = estac_file + '_n.txt'
        estac_file_z = estac_file + '_z.txt'

        # verifica que exista la estacion
        if ((os.path.isfile(estac_file_e) is False)):
            print("Faltan archivos de {}".format(estacion))
            continue

        # cargar modelo
        tiempo = np.loadtxt(estac_file_e, dtype=float, usecols=[0])
        desp_e = np.loadtxt(estac_file_e, dtype=float,  usecols=[1])
        desp_n = np.loadtxt(estac_file_n, dtype=float,  usecols=[1])
        desp_z = np.loadtxt(estac_file_z, dtype=float,  usecols=[1])
        estacion = [estaciones[i]]*len(tiempo)

        # nota: python pasa los valores a string en esta operacion
        # corregir usando .astype(np.float).T
        data_temp = np.array([estacion, tiempo,
                              desp_e, desp_n, desp_z],
                             dtype='str').T
        # limpiar datos
        data_temp2 = extraer_tiempo(data_temp, intervalo)

        # si se tienen menos de 10 mediciones descarta la estacion
        if len(data_temp2) < 10:
            continue

        # array con el nombre de la estacion de longitud nº estaciones
        lista_estac.append([estaciones[i], lon_est[i], lat_est[i]])

        # la primera vez data = data_temp, luego va apilando
        if data_vacio:
            data = data_temp2
            data_vacio = False
        else:
            data = np.vstack((data, data_temp2))

    # ###################################################################
    # OUTPUT
    # ##################################################################
    lista_final = np.array(lista_estac)
    return data, lista_final
