#! /usr/bin/env python3

import math
import numpy as np
import pdb
import os

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

import geometry as geo
import field as field
from example_figures import tensor_figure


# cargar malla
dir_path = os.path.dirname(os.path.realpath(__file__))
path_node = "{}/../example_files/malla.node".format(dir_path)
path_ele = "{}/../example_files/malla.ele".format(dir_path)
x_triangulos = np.loadtxt(path_node, usecols=[1], skiprows=1)
y_triangulos = np.loadtxt(path_node, usecols=[2], skiprows=1)
origen = [-76, -46]

# cargar los datos
dataset = '{}/../example_files/vect2007_2010'.format(dir_path)
lon_gps = np.loadtxt(dataset, usecols=[1])
lat_gps = np.loadtxt(dataset, usecols=[2])
v_n_gps = np.loadtxt(dataset, usecols=[4])
v_e_gps = np.loadtxt(dataset, usecols=[3])

# conversion mm/año a m/año
v_e_gps = v_e_gps/1000
v_n_gps = v_n_gps/1000

# proyectar malla y puntos gps
malla_proy = geo.geo2proj(x_triangulos, y_triangulos, origen[0], origen[1])
gps_proy = geo.geo2proj(lon_gps, lat_gps, origen[0], origen[1])

# interpolar en puntos de malla triangular
vy_tri = griddata(gps_proy, v_n_gps, malla_proy, method='cubic')
vx_tri = griddata(gps_proy, v_e_gps, malla_proy, method='cubic')

"""
NOTA FAKO
CON MALLA_PROY Y VY_TRI/VX_TRI PUEDES PLOTEAR LOS VECTORES
INTERPOLADOS EN LA MALLA
"""

# ############################################################################
# Funciones field()
# ############################################################################
gradiente = field.triangular_gradient(vx_tri, vy_tri, malla_proy[0],
                                      malla_proy[1], path_ele)
evalue, evector = field.principal_stress(gradiente)
# ############################################################################
# Graficar
# ############################################################################
fig = tensor_figure(x_triangulos, y_triangulos, evalue, evector)
# plt.show()
