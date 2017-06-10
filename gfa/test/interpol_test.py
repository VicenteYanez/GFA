#! /usr/bin/python3

import math
import numpy as np
import pdb
import os

import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import cartopy.crs as ccrs

from gfa.field_analysis import geometry as geo
from gfa.field_analysis import field
from gfa.field_analysis import FallaOkada as oka85
from gfa.field_analysis.fun_fig import vorticity_figure, cinematic_figure
from gfa.load_param import Config


def main():
    # make grid
    gridx, gridy = np.meshgrid(np.linspace(-74, -70.5, 50),
                               np.linspace(-46, -35, 100))
    gridx = gridx.ravel()
    gridy = gridy.ravel()

    # cargar los datos
    dir_path = os.path.dirname(os.path.realpath(__file__))
    dataset = '{}/../example_files/vect2007_2010'.format(dir_path)
    lon_gps = np.loadtxt(dataset, usecols=[1])
    lat_gps = np.loadtxt(dataset, usecols=[2])
    vn_gps = np.loadtxt(dataset, usecols=[4])
    ve_gps = np.loadtxt(dataset, usecols=[3])

    # conversion mm/año a m/año
    ve_gps = ve_gps/1000
    vn_gps = vn_gps/1000

    # cardozo&allmendinger strain
    b = np.reshape([ve_gps, vn_gps], (2*len(ve_gps), 1), order='F')
    gradiente = field.distance_weigthed2d(b, lon_gps, lat_gps, gridx, gridy,
                                          alfa=25000, dmin=100000)
    # tensor calculations
    S, W = field.velocitytensor_2d(gradiente[0], gradiente[1],
                                   gradiente[2], gradiente[3])
    wz = field.vertical_vorticity2d(W)
    Wk = field.cinematic_vorticity(S, W)

    # Graficar
    fig1 = vorticity_figure(gridx, gridy, wz)
    plt.plot(lon_gps, lat_gps, 'ko', transform=ccrs.PlateCarree())
    fname = "{}/../example_files/vorticity_weight_stat.png".format(dir_path)
    plt.savefig(fname, dpi=500, edgecolor='k', orientation='portrait',
                bbox_inches=None, pad_inches=0.1)
    fig2 = cinematic_figure(gridx, gridy, Wk)
    plt.plot(lon_gps, lat_gps, 'ko', transform=ccrs.PlateCarree())
    fname = "{}/../example_files/vorticityK_weight_stat.png".format(dir_path)
    plt.savefig(fname, dpi=500, edgecolor='k', orientation='portrait',
                bbox_inches=None, pad_inches=0.1)
    print('end function')


def okada_test():
    # geometria de falla esperable.
    dip = np.radians(90)    # [grados]
    W = 30000  # 197000.
    D = 10000+W*np.sin(dip)  # 60000.            # [metros] d en okada
    rake = np.radians(180)
    largo = 1100000.0
    strike = np.radians(10)

    # make grid
    gridx, gridy = np.meshgrid(np.linspace(-74, -70.5, 50),
                               np.linspace(-46, -35, 100))
    gridx = gridx.ravel()  # flattern the array
    gridy = gridy.ravel()

    # position stations
    dir_path = os.path.dirname(os.path.realpath(__file__))
    dataset = '{}/../example_files/vect2007_2010'.format(dir_path)
    lon_gps = np.loadtxt(dataset, usecols=[1])
    lat_gps = np.loadtxt(dataset, usecols=[2])

    # coordenadas geograficas del origen de okada
    latf = -48.0
    lonf = -74.0

    parametros = [D, W, largo, rake, strike, dip, latf, lonf]

    falla1 = oka85.FallaOkada(parametros, [lon_gps, lat_gps])
    Ue, Un, uz = falla1.desplaz()
    falla2 = oka85.FallaOkada(parametros, [gridx, gridy])
    oka_gradients = falla2.tensor_deformacion()

    # cardozo&allmendinger strain
    b = np.reshape([Ue, Un], (2*len(Ue), 1), order='F')
    gradiente = field.distance_weigthed2d(b, lon_gps, lat_gps, gridx, gridy,
                                          alfa=20000, dmin=100000)
    # tensor calculations
    S, W = field.velocitytensor_2d(gradiente[0], gradiente[1],
                                   gradiente[2], gradiente[3])
    wz = field.vertical_vorticity2d(W)
    Wk = field.cinematic_vorticity(S, W)

    # okada tensors calculation
    S_oka, W_oka = field.velocitytensor_2d(oka_gradients[0], oka_gradients[1],
                                           oka_gradients[2], oka_gradients[3])
    wz_oka = field.vertical_vorticity2d(W_oka)
    Wk_oka = field.cinematic_vorticity(S_oka, W_oka)

    # Graficar
    fig1 = vorticity_figure(gridx, gridy, wz)
    plt.plot(lon_gps, lat_gps, 'ko', transform=ccrs.PlateCarree())
    plt.quiver(lon_gps, lat_gps, Ue, Un)
    fname = "{}/../example_files/vorticity_weight_oka.png".format(dir_path)
    plt.savefig(fname, dpi=500, edgecolor='k', orientation='portrait',
                bbox_inches=None, pad_inches=0.1)
    fig2 = cinematic_figure(gridx, gridy, Wk)
    plt.plot(lon_gps, lat_gps, 'ko', transform=ccrs.PlateCarree())
    fname = "{}/../example_files/vorticityK_weight_oka.png".format(dir_path)
    plt.savefig(fname, dpi=500, edgecolor='k', orientation='portrait',
                bbox_inches=None, pad_inches=0.1)

    # figure with oka gradients
    fig3 = vorticity_figure(gridx, gridy, wz_oka)
    plt.plot(lon_gps, lat_gps, 'ko', transform=ccrs.PlateCarree())
    fname = "{}/../example_files/vorticity_okagradients.png".format(dir_path)
    plt.savefig(fname, dpi=500, edgecolor='k', orientation='portrait',
                bbox_inches=None, pad_inches=0.1)
    fig4 = cinematic_figure(gridx, gridy, Wk_oka)
    plt.plot(lon_gps, lat_gps, 'ko', transform=ccrs.PlateCarree())
    fname = "{}/../example_files/vorticityK_oka_gradients.png".format(dir_path)
    plt.savefig(fname, dpi=500, edgecolor='k', orientation='portrait',
                bbox_inches=None, pad_inches=0.1)
    print('end function')


main()
okada_test()
