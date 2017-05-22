#! /usr/bin/python3

"""
@author: Vicente Yáñez
"""

import os
import copy
import zipfile
import datetime as dt
import pdb

import numpy as np

from gfa.gnss_analysis.loadGPS import format_csn, format_model
import gfa.gnss_analysis.fun_vector


class TimeSeriesControl():
    """
    Clase con los metodos necesarios para cargar y graficar series de tiempo.

    Pensada para trabajar con bases de datos de muchas series de tiempo

    Para iniciar necesita la un archivo txt con la ubicacion
    de las estaciones, el directorio con las series de tiempo y
    un directorio donde se van a guardar los datos

    Variables globales
        Usuario
        tiempo

    Metodos:
        generar_modelo()
        pedir_estaciones()
        descargar_datos()

    """
    def __init__(self, codigo, lista_estac, dir_series, save_dir):
        self.dir_series = dir_series
        self.lista_estac = lista_estac
        self.savedir = '{}{}/'.format(save_dir, codigo)
        self.codigo = codigo
        self.clas = 'series'

        os.mkdir(self.savedir)

    def load_estations(self, lon_min, lon_max, lat_min, lat_max, tmin, tmax):
        """
        Metodo que recopila la información necesaria de las estaciones
        pedidas
        1. Carga el area y el intervalo de tiempo pedido por el usuario.
        2. Pide la selección de las estaciones

        Input
        lon_min, lon_max    : Rango de longitud solicitado
        lat_min, lat_max    : Rango de latitud solicitado
        tmin, tmax          : Rango de tiempo pedido

        Output
        nombre    : array-like con nombre de estaciones
        intervalo : intervalo de tiempo que contiene la estación
        error     : error de ajuste de la solucion actual
        """

        lon = [lon_min, lon_max]
        lat = [lat_min, lat_max]
        intervalo = [tmin, tmax]

        # Carga de datos
        data, lista = format_csn(lon, lat, intervalo, self.lista_estac,
                                 self.dir_series)

        self.data = data
        self.lista = lista

        return data, lista

    def savedata(self):
        """
        Guarda los datos cargados por load_estation
        """
        lista_estac = np.array(self.lista).T
        # guardar lista de estaciones
        head = 'estation,    longitude,    latitude'
        np.savetxt('{}{}_lista.txt'.format(self.savedir, self.clas),
                   self.lista, fmt='%s', header=head)

        # Crear directorio para almacenar los resultados
        save_series = '{}{}/'.format(self.savedir, self.clas)
        os.mkdir(save_series)

        # ciclos para extraer y guardar datos de cada estacion
        head = 'time[yr], Desp[mm] Est, North, Vertical, Err[mm] Est, North,\
Vertical'
        for estacion in self.data:
            save = np.array(estacion[1]).T
            np.savetxt('{}{}.txt'.format(save_series, estacion[0]), save,
                       fmt='%s',
                       header=head)
        return
