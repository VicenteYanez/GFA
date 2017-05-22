#! /usr/bin/python3

"""
@author: Vicente Yáñez
"""

import os
import copy
import json
import pdb

import numpy as np

from gfa.gnss_analysis.loadGPS import format_csn, format_model
from gfa.gnss_analysis.ModeloTrayectoria import ModeloTrayectoria
from gfa.gnss_analysis.TimeSeriesControl import TimeSeriesControl


class ModelControl(TimeSeriesControl):
    """
    Clase con los metodos necesarios para controlar
    las llamadas realizadas por una interfaz de
    usuario.
    Hereda metodos para manejar series de tiempo de TimeSeriesControl

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
        self.savedir = '{}{}/'.format(save_dir, codigo)
        self.codigo = codigo
        self.clas = 'modelo'
        self.lista_estac = lista_estac
        self.save_model = '{}Modelo/'.format(self.savedir)

        if os.path.isdir(self.savedir) is False:
            os.mkdir(self.savedir)
            os.mkdir(self.save_model)
            os.mkdir(self.save_model + 'comp/')

    def build_model_all(self, npoly, fperiods, eq_file):
        """
        Calculate the trajectory model for all the stations.
        Input:
            lista   : first list returned by loadGPS function
            data    : second list returned by loadGPS function
            npoly      : Degree of polynomial used by the Trajectory Model
            fperdios: periods of Furier function used by the Trajectory Model
            eq_file : file with the configuration of jumps/earthquakes
        """
        data = self.data
        lista = self.lista
        estac_pros = []

        for i, estac in enumerate(data):
            # posicion de estacion
            lon_estac = lista[i][1]
            lat_estac = lista[i][2]

            # Preludio Bevis
            model = ModeloTrayectoria(estac[1][0])
            tjump, tlt, tsc = self._load_events(estac[0], eq_file)
            model.tjump = tjump
            model.tlt = tlt
            model.tsc = tsc
            model.n = npoly
            model.fperiods = np.array(fperiods)

            # calcula y guarda modelo de trayectoria
            mdesp_e, res_e = model.modelo_trayectoria(estac[1][1])
            model.save_components(self.save_model, estac[0], '_e')
            modelo_e = model.total

            mdesp_n, res_n = model.modelo_trayectoria(estac[1][2])
            model.save_components(self.save_model, estac[0], '_n')
            modelo_n = model.total

            mdesp_u, res_u = model.modelo_trayectoria(estac[1][3])
            model.save_components(self.save_model, estac[0], '_z')
            modelo_u = model.total

            # save txt model files
            savearray = [estac[1][0], modelo_e, modelo_n, modelo_u]
            savearray = np.array(savearray).T
            head = 'time[yr], Est[mm], North[mm], Vertical[mm]'
            file_estac = '{}{}.txt'.format(self.save_model, estac[0])
            np.savetxt(file_estac, savearray, '%12.8f', header=head)

            # guardar txt con parametros
            parametros = {'polinomio': model.n,
                          'saltos': model.tjump.tolist(),
                          'Periodos Fourier': model.fperiods.tolist(),
                          'Inicio log': model.tlt.tolist(),
                          'Escala curva log': model.tsc.tolist()}
            residual = {'Residual E': res_e.tolist(),
                        'Residual N': res_n.tolist(),
                        'Residual Z': res_u.tolist()}

            parametros = json.dumps(parametros)
            residual = json.dumps(residual)

            # almacena las estaciones usadas
            pos_estac2 = [estac[0], lon_estac, lat_estac, parametros,
                          residual]
            estac_pros.append(pos_estac2)

        # guardar txt con lista de estaciones solucionadas
        head = 'estation,    longitude    latitude    Parameters    Error'
        np.savetxt('{}resume.txt'.format(self.savedir), estac_pros,
                   delimiter='    ', fmt='%s', header=head)

        return

    def load_model(self, lon_min, lon_max, lat_min, lat_max, tmin, tmax):
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
        data, lista = format_model(lon, lat, intervalo, self.lista_estac,
                                   self.dir_series)

        self.data = data
        self.lista = lista

        return data, lista

    def _load_events(self, name_est, earthq_file):
        """
        Function that calculate the jumps, tlt and tsc of the trajectory
        model, based in the config file earq_file.txt
        """
        # posicion de estacion
        try:
            lon_estacion = [x[:][1] for x in self.lista if x[:][0] == name_est]
            lat_estacion = [x[:][2] for x in self.lista if x[:][0] == name_est]
            # debe haber solo un dato de lon y lat para una estacion
            if len(lon_estacion) > 1 or len(lat_estacion) > 1:
                raise(ValueError)
        except ValueError:
            print('Hay más de una ubicación para una estacion en la db')
        finally:
            lon_estacion = lon_estacion[0]
            lat_estacion = lat_estacion[0]

        # terremotos
        fecha_evento = np.loadtxt(earthq_file, usecols=[2], skiprows=1)
        # intervalo de area de efecto de los eventos
        int_lat = np.loadtxt(earthq_file, usecols=[0, 1], skiprows=1)

        tjump = []
        tlt = []
        tsc = []

        # incorpora datos de sismos y saltos del archivo eventos.txt
        for n, evento in enumerate(fecha_evento):
            # si estacion esta en area de efecto del archivo eventos.txt
            if lat_estacion < int_lat[n][0] and lat_estacion > int_lat[n][1]:
                print("Incorpora sismo {} {} {}".format(name_est,
                                                        lat_estacion,
                                                        fecha_evento[n]))
                tjump.append(fecha_evento[n])
                tlt.append(fecha_evento[n])
                tsc.append(2)
            else:
                print("Ignora {} {} {} {}".format(name_est, lat_estacion,
                                                  int_lat[n][0],
                                                  int_lat[n][1]))
        tjump = np.array(tjump)
        tlt = np.array(tlt)
        tsc = np.array(tsc)

        return tjump, tlt, tsc

    def savedata(self):
        """
        Guarda los datos cargados por load_estation
        """
        lista_estac = np.array(self.lista).T
        # guardar lista de estaciones
        head = 'estation    longitude    latitude    Parameters    Error'
        np.savetxt('{}{}_lista.txt'.format(self.savedir, self.clas),
                   self.lista, fmt='%s', header=head, delimiter='    ')

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
