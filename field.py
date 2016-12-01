#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
@author: vicente
@date: December 2016

Conjunto de funciones relacionadas al campo
de velocidad.

"""

import numpy as np
from scipy import interpolate


def triangular_gradient(ux, uy, x, y, ele):
    """"
    Funcion que calcula el gradiente de una malla
    construida con el software triangle.
    Requiere que archivo .ele y .node esten en la misma carpeta

    Input: nombre de archivo .ele (1d narray)
    Output: ux/dx, ux/dy, uy/dx, uy/dy
    """

    # cargar datos de archivo. ele
    path_series = ele
    elem = np.loadtxt(path_series, usecols=[0], skiprows=1)
    vert_1 = np.loadtxt(path_series, usecols=[1], skiprows=1)
    vert_2 = np.loadtxt(path_series, usecols=[2], skiprows=1)
    vert_3 = np.loadtxt(path_series, usecols=[3], skiprows=1)

    # array con triangul ode referencia
    P_ref = np.array([[-1, 1, 0],
                      [-1, 0, 1]])

    # arrays de 0 para guardar los gradientes
    ungradiente = np.zeros(len(x))
    omega = np.zeros(len(x))
    gradiente = np.array([ungradiente, ungradiente,
                          ungradiente, ungradiente])

    # ciclo sobre los triangulos
    for i in range(len(elem)):
        # num nodos i-esimo triangulo
        nodos = [vert_1[i], vert_2[i], vert_3[i]]

        # componente de velocidad para los tres nodos del i-esimo triangulo
        uh_loc = [ux[nodos[0]-1], ux[nodos[1]-1],
                  ux[nodos[2]-1]]
        uh_loc += [uy[nodos[0]-1], uy[nodos[1]-1],
                   uy[nodos[2]-1]]

        # coordenadas de los nodos del i-esimo triangulo
        x = [x[nodos[0]-1], x[nodos[1]-1], x[nodos[2]-1]]
        y = [y[nodos[0]-1], y[nodos[1]-1], y[nodos[2]-1]]

        # matriz transpuesta e invertida B
        B = np.array([[x[1]-x[0], x[2]-x[0]],
                      [y[1]-y[0], y[2]-y[0]]])
        Bt = B.T
        Bt_inv = np.linalg.inv(Bt)
        gradP = np.dot(Bt_inv, P_ref)

        # calculo del area del i-esimo triangulo
        area = .5*(B[0][0]*B[1][1]-B[0][1]*B[1][0])

        # gradiente en el elemento
        DeP = np.zeros((4, 6))  # matriz misteriosa y magica
        DeP[:2, :3] = gradP
        DeP[2:, 3:] = gradP
        gradiente_elemento = np.dot(DeP, uh_loc)

        # sumatoria de gradientes para cada nodo
        for u in [0, 1, 2]:
            gradiente[0][nodos[u]-1] += gradiente_elemento[0]*area
            gradiente[1][nodos[u]-1] += gradiente_elemento[1]*area
            gradiente[2][nodos[u]-1] += gradiente_elemento[2]*area
            gradiente[3][nodos[u]-1] += gradiente_elemento[3]*area
            omega[nodos[u]-1] += area
    # division elemento a elemento entre gradientes y area
    gradiente[0] = np.divide(gradiente[0], omega)  # grad uxdx
    gradiente[1] = np.divide(gradiente[1], omega)  # grad uxdy
    gradiente[2] = np.divide(gradiente[2], omega)  # grad uydx
    gradiente[3] = np.divide(gradiente[3], omega)  # grad uydy

    return gradiente


def vorticidad_general(self, uxdx, uxdy, uydx, uydy):
    S2 = uxdx**2 + uydy**2
    S2 += 0.5*(uxdy + uydx)**2
    S2 = np.sqrt(S2)

    # vorticidad vertical
    wz = uydx - uxdy

    # invariant of vorticity
    W2 = np.sqrt(0.5*(wz**2))

    # kinematic vorticity.
    Wk = W2/S2
    Wk = np.ma.masked_array(Wk, mask=np.isnan(Wk))  # ignorar valores nan

    return S2, wz, Wk


def traza_S(self, uxdx, uydy):
    S = uxdx + uydy

    return S


def vorticidad_monoclinica(self, uxdy, uydy):
    """
    En vorticidad monoclinica eje x es
    paralelo a la falla. Y deformación en este eje
    se asume como 0.
    uxdx, uydx = 0
    """
    S2 = uydy**2
    S2 += 0.5*(uxdy)**2
    S2 = np.sqrt(S2)

    # vorticidad vertical
    wz = uxdy

    # invariant of vorticity
    W2 = np.sqrt(0.5*(wz**2))

    # kinematic vorticity.
    Wk = W2/S2
    Wk = np.ma.masked_array(Wk, mask=np.isnan(Wk))  # ignorar valores nan

    return S2, wz, Wk

def limpiar_Wk(self, S2, Wk):
    """
    Elimina valores de Wk que pueden haber sido elebados por un S
    muy bajo
    """
    for i in range(len(S2)):
        if S2[i] < 0.1:
            Wk[i] = float('nan')
        if Wk[i] > 5:
            Wk[i] = float('nan')
    return Wk
