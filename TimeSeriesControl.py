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

from loadGPS import format_csn
import graficar_serie
import fun_tools
import fun_vector


class TimeSeriesControl():
    """
    Clase con los metodos necesarios para cargar y graficar series de tiempo

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
        self.savedir = save_dir
        self.codigo = codigo

    def load_estations(self, lon_min, lon_max, lat_min, lat_max, tmin, tmax):
        """
        Metodo que recopila la información necesaria de las estaciones
        pedidas
        1. Carga el area y el intervalo de tiempo pedido por el usuario.
        2. Pide la selección de las estaciones
        3. Pide la info de la solucion a la DB.

        Input
        Codigo    : Identificador de consulta, puede ser cualquier string

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
        # Crear directorio para almacenar los resultados
        save_dir = '{}{}_series/'.format(self.savedir, self.codigo)
        os.mkdir(save_dir)

        # ciclos para extraer y guardar datos de cada estacion
        for estacion in self.data:
            save = np.array(estacion[1]).T
            np.savetxt('{}{}.txt'.format(save_dir, estacion[0]), save, fmt='%s')
        return

    def to_graph(self, codigo, estacion, *recta):
        """
        Llama a función para crear y un gŕafico de una serie de tiempo.
        Se puede indicar el nombre de una estacion en vez de tener
        que leerla desde peticion_graficos.txt

        Los graficos se construyen desde:
            Tabla estaciones /proj/tmp/tabla_estaciones.txt
            Series de tiempo /proj/tmp/series/
            Modelo /proj/tmp/modelo/

        Input
        Codigo fecha

        Output
        archivos .jpg guardado en /temp/gps_graficos

        """
        # ####################################################################
        # CARGA DATOS
        # ####################################################################
        file_serie = self.tmp + codigo + '_series/' + estacion + '.txt'
        file_modelo = self.tmp + codigo + '_modelo/' + estacion + '.txt'
        serie = np.loadtxt(file_serie, usecols=range(1, 5)).T
        modelo = np.loadtxt(file_modelo, usecols=range(1, 5)).T

        # ####################################################################
        # CONSTRUCCIÓN DEL GRÁFICO
        # ####################################################################
        output = self.tmp + codigo + '_gps_graficos/'
        if os.path.isdir(output) is False:
            os.mkdir(output)
        grafico = Grafico(output)

        # Incluye valores de recta si esta fue dada
        if recta:
            grafico.graficar(estacion, serie[0], serie[1:4], modelo[1:4],
                             False, recta)
        else:
            grafico.graficar(estacion, serie[0], serie[1:4], modelo[1:4],
                             False)

        return

    def pedir_descarga(self, codigo):
        """
        Agrupa los datos recopilados y trabajados para
        su descarga.
        Input:
            codigo
        """
        # Carpetas a guardar
        dir_serie = codigo + '_series/'
        dir_vectores = codigo + '_vectores/'
        dir_graficos = codigo + '_gps_graficos/'
        file_param = codigo + '_tabla_estaciones.txt'

        # crear zip
        zipf = zipfile.ZipFile(self.tmp + codigo + '.zip', 'w',
                               zipfile.ZIP_DEFLATED)

        self._zipdir(self.tmp + dir_serie, dir_serie, zipf)
        self._zipdir(self.tmp + dir_modelo, dir_modelo, zipf)
        self._zipdir(self.tmp + dir_vectores, dir_vectores, zipf)
        self._zipdir(self.tmp + dir_graficos, dir_graficos, zipf)
        self._zipdir(self.tmp + file_param, file_param, zipf)

        zipf.close()

        return

    def _zipdir(self, path, dir_, ziph):
        # ziph is zipfile handle
        for root, dirs, files in os.walk(path):
            for file_ in files:
                ziph.write(os.path.join(root, file_), dir_ + file_)
