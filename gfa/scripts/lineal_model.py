#! /usr/bin/env python3

"""
@author: vicente

Script para realizar ajustes lineales de series de tiempo de datos GPS
y permite al usuario descartar o modificar la serie de tiempo individualmente.

Argumentos:
0: archivo txt con nombre y coordenadas de estaciones
1: carpeta donde se ubican las series de tiempo
2: Tiempo Mínimo
3: Tiempo Máximo

Format estaciones.txt
Nombre | Lon | Lat | Ve | Vn | Cov e | Cov n

Formato Series de Tiempo
Tiempo  Vel Error

Nombre de archivos de series de tiempo:
CODIGOgps_Componente.out
Ej:
AGUIgps_e.out
"""

import os
import sys
import pdb
from math import fabs, floor, ceil

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


# funcion de desplazamiento en funcion del tiempo
def modelo(tiempo, velocidad, posicion_inicial):
    return tiempo * velocidad + posicion_inicial


def ajustar(est, tiempo_e, tiempo_n, desp_e, desp_n, desp_e_err, desp_n_err,
            guardar=False):
    """
    Funcion que sirve para plotear una estacion
    con su ajuste definido por la funcion modelo

    tiempo, desp y error de desplazamiento son arrays con la series de tiempo
    completa para la estación 'est'
    """
    # funcion para imprimir recta que se ajusta a la curva
    # valor_e[0] es la velocidad, y cov_e[0] error en la vel
    valor_e, cov_e = curve_fit(modelo, tiempo_e, desp_e)
    valor_n, cov_n = curve_fit(modelo, tiempo_n, desp_n)

    # calculo del modelo
    modelo_e = [modelo(t, valor_e[0], valor_e[1]) for t in tiempo_e]
    modelo_n = [modelo(t, valor_n[0], valor_n[1]) for t in tiempo_n]

    fig = plt.figure(figsize=(20, 8))
    plt.ticklabel_format(useOffset=False)
    plt.errorbar(tiempo_e, desp_e, yerr=desp_e_err, color='r',
                 linestyle='None', label='Este')
    plt.errorbar(tiempo_n, desp_n, yerr=desp_n_err, color='g',
                 linestyle='None', label='Norte')
    plt.plot(tiempo_e, desp_e, 'ro', tiempo_e, modelo_e, '-k')
    plt.plot(tiempo_n, desp_n, 'go', tiempo_n, modelo_n, '-k')
    # formato
    year_min = int(floor(np.array(tiempo_e).min()))
    year_max = int(ceil(np.array(tiempo_e).max()))
    xticks = range(year_min, year_max, 1)

    plt.xticks(xticks, fontsize=14)

    plt.xlabel('tiempo [años]', fontsize=12, weight="bold")
    plt.ylabel('desplazamiento [mm]', fontsize=14, weight="bold")
    # adornos
    titulo = "Estacion {0} | Velocidad Este {1:4.2f} mm/año | velocidad Norte: \
{2:4.2f} mm/año".format(est, valor_e[0], valor_n[0])
    plt.title(titulo, fontsize=13)
    plt.show(block=True)

    # salvar grafico
    if guardar is True:
        plt.savefig(est)

    return valor_e, cov_e, valor_n, cov_n


def main(path_estaciones, path_series):
    tmin = float(sys.argv[3])
    tmax = float(sys.argv[4])

    # Carga de ubicacion de estaciones
    estaciones = np.loadtxt(path_estaciones, usecols=[0], dtype='S5')
    lon = np.loadtxt(path_estaciones, usecols=[2])
    lat = np.loadtxt(path_estaciones, usecols=[1])

    # array para guardar las estaciones con sus velocidades
    vel_estaciones = []
    estaciones_ausentes = []
    for i, estacion in enumerate(estaciones):
        # revisar si el archivo existe y si no es vacio
        estacion = estacion.decode('utf-8')
        e_file = "{}{}gps_e.out".format(path_series, estacion)
        n_file = "{}{}gps_n.out".format(path_series, estacion)
        e_exist = os.path.isfile(e_file)
        n_exist = os.path.isfile(n_file)

        if e_exist and n_exist:
            stinfo = os.stat(e_file)
            if stinfo.st_size > 60:
                # abrir archivos de estacion y carga de datos
                data_e = np.loadtxt(e_file)
                data_n = np.loadtxt(n_file)
                tiempo_e = data_e.T[0]
                tiempo_n = data_n.T[0]
                desp_e = data_e.T[1]
                desp_n = data_n.T[1]
                desp_e_err = data_e.T[2]
                desp_n_err = data_n.T[2]

                # plotear solo rangos de tiempo
                T_e = []
                T_n = []
                D_e = []
                D_n = []
                Derr_e = []
                Derr_n = []

                for u, valor in enumerate(tiempo_e):
                    if tiempo_e[u] >= tmin and tiempo_e[u] <= tmax:
                        T_e.append(tiempo_e[u])
                        T_n.append(tiempo_n[u])
                        D_e.append(desp_e[u])
                        D_n.append(desp_n[u])
                        Derr_e.append(desp_e_err[u])
                        Derr_n.append(desp_n_err[u])
                if len(T_e) < 2:
                    print("No hay datos para {}".format(estacion))
                    a = "No hay datos para {}".format(estacion)
                    estaciones_ausentes.append(a)
                else:
                    # plot de movimiento de estacion con su rango de error
                    valor_e, cov_e, valor_n, cov_n = ajustar(estacion,
                                                             T_e, T_n, D_e, D_n,
                                                             Derr_e, Derr_n)
                    ToF = input("¿Incluir estacion {}? Si, Editar, No (S, E, [N])\
    \n".format(estacion))
                    if ToF == "S":
                        # guardar valores
                        vel_estaciones.append(np.array([estacion, lon[i], lat[i],
                                                        valor_e[0], valor_n[0],
                                                        cov_e[0][0], cov_n[0][0]]))
                    elif ToF == "E":
                        # eliminar datos con errores altos
                        errlim = float(input("Error no permitido \n"))
                        T_e2 = []
                        T_n2 = []
                        D_e2 = []
                        D_n2 = []
                        Derr_e2 = []
                        Derr_n2 = []
                        for u in range(len(T_e)):
                            if Derr_e[u] < errlim:
                                T_e2.append(tiempo_e[u])
                                D_e2.append(desp_e[u])
                                Derr_e2.append(desp_e_err[u])
                            if Derr_n[u] < errlim:
                                T_n2.append(tiempo_n[u])
                                D_n2.append(desp_n[u])
                                Derr_n2.append(desp_n_err[u])
                        # recalcular y plotear grafico estación
                        if len(T_e2) > 2:
                            valor_e, cov_e, valor_n, cov_n = ajustar(estacion,
                                                                     T_e2, T_n2,
                                                                     D_e2, D_n2,
                                                                     Derr_e2,
                                                                     Derr_n2)
                        else:
                            print("Hay más error que: ", errlim)
                        # elimar un rango de tiempo poco preciso
                        cont = input("¿Eliminar rango de tiempo?.S/N \n")
                        if cont == "S":
                            fecha_min = float(input("Eliminar desde: \n"))
                            fecha_max = float(input("Eliminar hasta: \n"))
                            T_e3 = []
                            T_n3 = []
                            D_e3 = []
                            D_n3 = []
                            Derr_e3 = []
                            Derr_n3 = []
                            # adjuntar lo que no se encuentre en el rango de tiempo
                            for u in range(len(T_e2)):
                                if T_e2[u] < fecha_min or T_e2[u] > fecha_max:
                                    T_e3.append(T_e2[u])
                                    D_e3.append(D_e2[u])
                                    Derr_e3.append(Derr_e2[u])
                                if T_n2[u] < fecha_min or T_e2[u] > fecha_max:
                                    T_n3.append(T_n2[u])
                                    D_n3.append(D_n2[u])
                                    Derr_n3.append(Derr_n2[u])
                            # recalcular y plotear grafico estación
                            if len(T_e3) > 2:
                                valor_e, cov_e, valor_n, cov_n = ajustar(estacion,
                                                                         T_e3,
                                                                         T_n3,
                                                                         D_e3,
                                                                         D_n3,
                                                                         Derr_e3,
                                                                         Derr_n3)
                        ToF = input("¿Incluir estacion? Si, No (S, [N])\n \
                                    ".format(estacion))
                        if ToF == "S":
                            vel_estaciones.append(np.array([estacion,
                                                           lon[i], lat[i],
                                                           valor_e[0], valor_n[0],
                                                           cov_e[0][0],
                                                           cov_n[0][0]]))
                        else:
                            print("Estación no incorporada ", estacion)
                            estaciones_ausentes.append(np.array([estacion,
                                                                lon[i], lat[i],
                                                                valor_e[0],
                                                                valor_n[0],
                                                                cov_e[0][0],
                                                                cov_n[0][0]]))
                    else:
                        print("Estación no incorporada ", estacion)
                        estaciones_ausentes.append(np.array([estacion,
                                                            lon[i], lat[i],
                                                            valor_e[0], valor_n[0],
                                                            cov_e[0][0],
                                                            cov_n[0][0]]))
            else:
                print("Archivo vacio ", estacion)
                a = "Archivo vacio {}".format(estacion)
                estaciones_ausentes.append(a)
        else:
            print("No existe archivo de ", estacion)
            a = "No existe archivo de {}".format(estacion)
            estaciones_ausentes.append(a)

    np.savetxt("vel_estaciones", vel_estaciones, fmt="%s")
    np.savetxt("estaciones_ausentes", estaciones_ausentes, fmt="%s")
    print("Fin")
