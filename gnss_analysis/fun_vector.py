#! /usr/bin/python3


import numpy as np
from scipy.optimize import curve_fit

"""
Colección de funciones para calcular un vector desde
el modelo de trayectoria
"""


def fun_lineal(x, b, c):
    return x * b + c


def tangente(ti, t, e, n, z):
    """
    Genera un vector en el tiempo ti
    para un modelo de desplazamiento (n, e, z)
    con un tiempo t

    Input:
    ti: ubicación de la recta tangente
    t:  array con tiempo del modelo
    e:  array con desplazamiento este
    n:  array con desplazamiento norte
    z:  array con desplazamiento vertical
    """
    # evalua la velocidad en base a los dos puntos mas cercanos del modelo
    idx = (np.abs(t-ti)).argmin()

    # dependiendo de la posicion de idx, considera el punto mas cercano
    # a la izquierda o a la derecha
    if idx-1 > 0:
        vel_e = (e[idx] - e[idx-1]) / (t[idx]-t[idx-1])
        vel_n = (n[idx] - n[idx-1]) / (t[idx]-t[idx-1])
        vel_z = (z[idx] - z[idx-1]) / (t[idx]-t[idx-1])
        # y = mx + c
        # calculo de c
        c_e = vel_e*(-t[idx-1]) + e[idx-1]
        c_n = vel_n*(-t[idx-1]) + n[idx-1]
        c_z = vel_z*(-t[idx-1]) + z[idx-1]

    else:
        vel_e = (e[idx+1] - e[idx]) / (t[idx+1]-t[idx])
        vel_n = (n[idx+1] - n[idx]) / (t[idx+1]-t[idx])
        vel_z = (z[idx+1] - z[idx]) / (t[idx+1]-t[idx])
        # y = mx + c
        # calculo de c
        c_e = vel_e*(-t[idx]) + e[idx]
        c_n = vel_n*(-t[idx]) + n[idx]
        c_z = vel_z*(-t[idx]) + z[idx]

    return np.array([vel_e, vel_n, vel_z]), np.array([c_e, c_n, c_z])


def fit(intervalo, t, e, n, z):
    """
    Genera un vector realizando un ajuste lineal
    sobre una serie en un intervalo de tiempo
    """
    # obtiene datos de intervalo de tiempo
    id_t = (t >= intervalo[0]) & (t <= intervalo[1])
    t = t[id_t]
    e = e[id_t]
    n = n[id_t]
    z = z[id_t]
    if len(t) < 10:
        print("Estacion sin datos para calcular una vel media")
        return False

    # ajuste de datos
    pva_e, cov_e = curve_fit(fun_lineal, t, e)
    pva_n, cov_n = curve_fit(fun_lineal, t, n)
    pva_z, cov_z = curve_fit(fun_lineal, t, z)

    err_e = np.sqrt(np.diag(cov_e))
    err_n = np.sqrt(np.diag(cov_n))
    err_z = np.sqrt(np.diag(cov_z))

    return (np.array([pva_e[0], pva_n[0], pva_z[0]]),
            np.array([pva_e[1], pva_n[1], pva_z[1]]),
            np.array([err_e[0], err_n[0], err_z[0]]))
