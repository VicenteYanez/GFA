#! /usr/bin/env python3

"Function to plot the examples script"

import pdb

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata


def greater_absolute(value1, value2):
    if abs(value1) >= abs(value2):
        return np.array([-value1, value1])
    else:
        return np.array([value2, -value2])


def vorticity_figure(x, y, z):
    """
    Figure of example_vorticity.py
    """
    fig = plt.figure(figsize=(6, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())

    paralelos = np.arange(-46, -35, 2)
    meridianos = np.arange(-76, -70, 2)
    labelE = [1, 0, 0, 0]
    labelS = [1, 0, 0, 1]
    grosor_linea = 0.5

    # extension of the map
    latmin_m = y.min() - 0.5
    latmax_m = y.max() + 0.5
    lonmin_m = x.min() - 0.5
    lonmax_m = x.max() + 0.5
    #  zone to interpol
    latmin = y.min()
    latmax = y.max()
    lonmin = x.min()
    lonmax = x.max()

    ax.set_extent([lonmin_m, lonmax_m, latmin_m, latmax_m])
    ax.set_xticks(meridianos, crs=ccrs.PlateCarree())
    ax.set_yticks(paralelos, crs=ccrs.PlateCarree())
    # ax.gridlines(crs=ccrs.PlateCarree())

    x2 = []
    y2 = []
    z2 = []
    for i, value in enumerate(z):
        if np.isnan(z[i]) is False:
            x2.append(x[i])
            y2.append(y[i])
            z2.append(z[i])

    # interpolation grid
    resolucion = 400
    lats = np.linspace(latmin, latmax, resolucion)
    lons = np.linspace(lonmin, lonmax, resolucion)
    y_mapa, x_mapa = np.meshgrid(lats, lons)
    z_mapa = griddata((x, y), z, (x_mapa, y_mapa), method='linear')

    # mask nan values
    z_mapa = np.ma.masked_array(z_mapa, mask=np.isnan(z_mapa))

    # plot grid
    cmap = plt.cm.seismic
    cmap.set_under(color="black", alpha="0")
    # colorbar scale
    z_mapa_esc = z_mapa*360  # to convert to degrees/year
    vesc = greater_absolute(z_mapa_esc.min(), z_mapa_esc.max())
    im = ax.pcolormesh(x_mapa, y_mapa, z_mapa_esc,
                       cmap=cmap,
                       transform=ccrs.PlateCarree(),
                       vmin=vesc.min(), vmax=vesc.max(),
                       alpha=1, zorder=1)
    cb = plt.colorbar(im)
    cb.set_label('º/year')

    # final details
    ax.coastlines(resolution='10m')
    return fig


def cinematic_figure(x, y, z):
    fig = plt.figure(figsize=(6, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())

    paralelos = np.arange(-46, -35, 2)
    meridianos = np.arange(-76, -70, 2)
    labelE = [1, 0, 0, 0]
    labelS = [1, 0, 0, 1]
    grosor_linea = 0.5

    # extension of the map
    latmin_m = y.min() - 0.5
    latmax_m = y.max() + 0.5
    lonmin_m = x.min() - 0.5
    lonmax_m = x.max() + 0.5
    #  zone to interpol
    latmin = y.min()
    latmax = y.max()
    lonmin = x.min()
    lonmax = x.max()

    ax.set_extent([lonmin_m, lonmax_m, latmin_m, latmax_m])
    ax.set_xticks(meridianos, crs=ccrs.PlateCarree())
    ax.set_yticks(paralelos, crs=ccrs.PlateCarree())
    # ax.gridlines(crs=ccrs.PlateCarree())

    x2 = []
    y2 = []
    z2 = []
    for i, value in enumerate(z):
        if np.isnan(z[i]) is False:
            x2.append(x[i])
            y2.append(y[i])
            z2.append(z[i])

    # interpolation grid
    resolucion = 1000
    lats = np.linspace(latmin, latmax, resolucion)
    lons = np.linspace(lonmin, lonmax, resolucion)
    y_mapa, x_mapa = np.meshgrid(lats, lons)
    z_mapa = griddata((x, y), z, (x_mapa, y_mapa), method='linear')

    # mask nan values
    z_mapa = np.ma.masked_array(z_mapa, mask=np.isnan(z_mapa))

    # plot grid
    cmap = plt.cm.jet
    cmap.set_under(color="black", alpha="0")

    im = ax.pcolormesh(x_mapa, y_mapa, z_mapa,
                       cmap=cmap,
                       transform=ccrs.PlateCarree(),
                       vmin=0, vmax=1,
                       alpha=1, zorder=1)
    cb = plt.colorbar(im)
    cb.set_label('Cinematic Vorticity')

    # final details
    ax.coastlines(resolution='10m')
    return fig


def add_vector(figure, axes, x, y, vx, vy):
    figure.plot(x, y, vx, vy)

    return figure, axes



def tensor_figure(x, y, evalue, evector):
    """
    Figure of example_principalstress.py
    """
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())

    paralelos = np.arange(-46, -35, 2)
    meridianos = np.arange(-76, -70, 2)
    labelE = [1, 0, 0, 0]
    labelS = [1, 0, 0, 1]
    grosor_linea = 0.5

    # extension of the map
    latmin_m = y.min() - 0.5
    latmax_m = y.max() + 0.5
    lonmin_m = x.min() - 0.5
    lonmax_m = x.max() + 0.5
    #  zone to interpole
    latmin = y.min()
    latmax = y.max()
    lonmin = x.min()
    lonmax = x.max()

    ax.set_extent([lonmin_m, lonmax_m, latmin_m, latmax_m])
    ax.set_xticks(meridianos, crs=ccrs.PlateCarree())
    ax.set_yticks(paralelos, crs=ccrs.PlateCarree())
    ax.gridlines(crs=ccrs.PlateCarree())

    x2 = []
    y2 = []
    for i, value in enumerate(x):
        if np.isnan(x[i]) is False:
            x2.append(x[i])
            y2.append(y[i])

    # eigen-vector and eigen-value format
    xf = [i for i in np.array(evector)[:, 0, 0] if i.imag != 0]
    yf = np.array(evector)[:, 0, 1]
    print(xf)
    pdb.set_trace()
    # plot tensor
    plt.quiver(x, y, evector.T[0], color='red')  # first vector
    plt.quiver(x, y, color='blue')  # second vector
    plt.quiver(x, y, evector.T[0], color='red')  # oppose first vector
    plt.quiver(x, y, color='blue')  # oppose second vector

    # final details
    ax.coastlines(resolution='10m')
    return fig
