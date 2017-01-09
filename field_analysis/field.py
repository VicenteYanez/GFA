#! /usr/bin/python3

"""
@author: vicente
@date: December 2016

Set of functions related with the velocity field
"""

import numpy as np


def triangular_gradient(ux, uy, x, y, ele):
    """"
    Function that calculate the gradient in a mesh, build with
    triangle software
    Need that the .ele and .node files are in the same directory.

    Input:
        ux, uy : velocity in the x and y axis
        x, y   : triangular mesh
        ele: path of the .ele file
    Output:
        ux/dx, ux/dy, uy/dx, uy/dy
    """

    # cargar datos de archivo. ele
    path_series = ele
    elem = np.loadtxt(path_series, usecols=[0], skiprows=1, dtype=int)
    vert_1 = np.loadtxt(path_series, usecols=[1], skiprows=1, dtype=int)
    vert_2 = np.loadtxt(path_series, usecols=[2], skiprows=1, dtype=int)
    vert_3 = np.loadtxt(path_series, usecols=[3], skiprows=1, dtype=int)
    # array con triangulo de referencia
    P_ref = np.array([[-1, 1, 0],
                      [-1, 0, 1]])

    # arrays de 0 para guardar los gradientes
    ungradiente = np.zeros(len(x))
    omega = np.zeros(len(x))
    gradiente = np.array([ungradiente, ungradiente,
                          ungradiente, ungradiente])

    # ciclo sobre los triangulos
    for i, elemento in enumerate(elem):
        # num nodos i-esimo triangulo
        nodos = [vert_1[i], vert_2[i], vert_3[i]]

        # componente de velocidad para los tres nodos del i-esimo triangulo
        uh_loc = [ux[nodos[0]-1], ux[nodos[1]-1], ux[nodos[2]-1]]
        uh_loc += [uy[nodos[0]-1], uy[nodos[1]-1], uy[nodos[2]-1]]

        # coordenadas de los nodos del i-esimo triangulo
        xi = [x[nodos[0]-1], x[nodos[1]-1], x[nodos[2]-1]]
        yi = [y[nodos[0]-1], y[nodos[1]-1], y[nodos[2]-1]]

        # matriz transpuesta e invertida B
        B = np.array([[xi[1]-xi[0], xi[2]-xi[0]],
                      [yi[1]-yi[0], yi[2]-yi[0]]])
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


def velocitytensor_2d(uxdx, uxdy, uydx, uydy):
    """
    Function that calculates the elongation and vorticity tensors,
    from the velocity gradients for two dimensions.

    The returning tensors have this form

    T = T11 T12
        T21 T22
    where Txy correspond to T[x][y]
    """
    # convierte input a array
    uxdx = np.array(uxdx)
    uxdy = np.array(uxdy)
    uydx = np.array(uydx)
    uydy = np.array(uydy)

    tensorS = 0.5*np.array([[2*uxdx, uxdy+uydx],
                            [uydx+uxdy, 2*uydy]])
    tensorW = 0.5*np.array([[uxdx-uxdx, uxdy-uydx],
                            [uydx-uxdy, uydy-uydy]])

    return tensorS, tensorW


def velocitytensor_3d(uxdx, uxdy, uxdz, uydx, uydy, uydz, uzdx, uzdy, uzdz):
    """
    Function that calculates the elongation and vorticity tensors,
    from the velocity gradients for tree dimensions.

    The returning tensors have this form

    S = S11 S12 S13
        S21 S22 S23
        S31 S32 S33
    where Txy correspond to T[x][y]
    """
    # convierte input a array
    uxdx = np.array(uxdx)
    uxdy = np.array(uxdy)
    uxdz = np.array(uxdz)
    uydx = np.array(uydx)
    uydy = np.array(uydy)
    uydz = np.array(uydz)
    uzdx = np.array(uzdx)
    uzdy = np.array(uzdy)
    uzdz = np.array(uzdz)

    tensorS = 0.5*np.array([[2*uxdx, uxdy+uydx, uxdz+uzdx],
                            [uydx+uxdy, 2*uydy, uydz+uzdy],
                            [uzdx+uxdz, uzdy+uydz, 2*uzdz]])

    tensorW = 0.5*np.array([[uxdx-uxdx, uxdy-uydx, uxdz-uzdx],
                            [uydx-uxdy, uydy-uydy, uydz+uzdy],
                            [uzdx-uxdz, uzdy-uydz, uzdz-uzdz]])

    return tensorS, tensorW


def vertical_vorticity2d(tensorW):
    """
    Calculate the vorticity vector from the vorticity tensor in two dimension
    """
    wz = 2*tensorW[1][0]
    return wz


def vorticity_vector3d(tensorW):
    """
    Calculate the vorticity vector from the vorticity tensor in
    tree dimensions
    """
    w = 2*np.array([tensorW[2][1], tensorW[0][2], tensorW[1][0]])

    return w


def cinematic_vorticity(tensorS, tensorW):
    module_S = 0
    module_W = 0
    for s1, w1 in zip(tensorS, tensorW):
        for s2, w2 in zip(s1, w1):
            module_S += s2**2
            module_W += w2**2

    module_S = np.sqrt(module_S)
    module_W = np.sqrt(module_W)

    Wk = module_S/module_W
    Wk = np.ma.masked_array(Wk, mask=np.isnan(Wk))  # ignorar valores nan

    return Wk


def traza_tensor(tensor, grado):
    traza = 0
    for i, valor in enumerate(tensor):
        traza += valor[i]**grado
    return traza


def principal_stress(tensor):
    return None


def clean_array(array_pp, value_min, value_max, param_=()):
    """
    Function that replace with nan, the values of array_pp that
    are outside the interval value_min - value_max

    Optionally, you can give one o more list or numpy.arrays to
    remove the same indices of array_pp
    """
    # copia de arrays entregados para evitar modificacion de mutables
    copy_app = array_pp[:]
    len_pp = len(copy_app)
    param_ = param_[:]

    # formato de param_ a una lista
    if param_ and type(param_) is not tuple:
        param_ = [param_]
    elif type(param_) is tuple:
        param_ = list(param_)

    for i, value in enumerate(copy_app):
        if value < value_min or value > value_max:
            copy_app[i] = float('nan')

            # ciclo sobre los arrays opcionales
            if param_:
                for array in param_:
                    if len(array) != len_pp:
                        raise ValueError('Invalid inputs: input arrays has \
                                         not the same lenght')
                    else:
                        array[i] = float('nan')
    if param_:
        return copy_app, param_
    else:
        return copy_app
