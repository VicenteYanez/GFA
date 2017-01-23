#! /usr/bin/python3

"""
@author: Vicente Yáñez
"""

import os
import copy
import json
import zipfile
import datetime as dt
import pdb

import numpy as np

from loadGPS import format_csn, format_model
from ModeloTrayectoria import ModeloTrayectoria
from TimeSeriesControl import TimeSeriesControl
import graficar_serie
import fun_vector


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

    def build_model(self, data, param):
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
        np.savetxt('{}resumen_solucion.txt'.format(self.savedir), estac_pros,
                   delimiter='    ', fmt='%s', header=head)

        return

    def calc_vector(self, estacion, tipo_vector, *t0):
        """
        Metodo que llama a una funcion para crear un vector
        dependiendo del tipo de calculo elegido
        Input
        *t0: ubicacion de la recta tangente
        """
        # ####################################################################
        # CARGA DATOS
        # ####################################################################
        file_modelo = self.tmp + self.codigo + '_modelo/' + estacion + '.txt'
        modelo = np.loadtxt(file_modelo, usecols=(1, 2, 3, 4)).T

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

        # ####################################################################
        # GUARDAR VECTOR
        # ####################################################################
        # Crear directorio para almacenar vectores
        vector_dir = self.tmp + self.codigo + '_vectores/'
        vector_file = vector_dir + 'vectores'
        if os.path.isdir(vector_dir) is False:
            os.mkdir(vector_dir)

        # incorporar nombre de estacion a array para guardar
        save_vector = np.hstack(([estacion], vector))

        # guardar vector estacion
        if os.path.isfile(vector_file):
            vectores = np.loadtxt(vector_file, dtype=str)
            vectores = np.vstack((vectores, save_vector))
            np.savetxt(vector_file, vectores, fmt='%s', delimiter='    ')
        else:
            np.savetxt(vector_file, [save_vector], fmt='%s', delimiter='    ')

        # ####################################################################
        # RECONSTRUIR GRAFICO
        # ####################################################################
        self.pedir_grafico(self.codigo, estacion, vector, c)

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

    def _actualizar_resumen(self, estacion, parametros, residual,
                            autor='auto'):
        """
        Actualiza los valores en la tabla de estaciones generados por
        main_seleccionar
        """
        archivo = self.tmp + self.codigo+'_tabla_estaciones.txt'
        data = np.loadtxt(archivo, dtype=str, delimiter='    ')

        # Modificar
        for i in range(len(data)):
            if data[i][0] == estacion:
                data[i] = np.array([estacion, data[i][1], data[i][2],
                                    parametros, residual])

        # actualizar txt
        np.savetxt(archivo, data, fmt='%s', delimiter='    ')
        return

    def _load_events(self, name_est, earthq_file):

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
