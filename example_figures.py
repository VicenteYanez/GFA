#! /usr/bin/env python3

"Funciones que crean las figuras de los ejemplos"


import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata


def vorticity_figure(x, y, z):
    """
    Figura de example.vorticity.py
    """
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())

    # preludio
    paralelos = np.arange(-46, -35, 2)
    meridianos = np.arange(-76, -70, 2)
    labelE = [1, 0, 0, 0]
    labelS = [1, 0, 0, 1]
    grosor_linea = 0.5

    # rango de coordenadas del mapa
    latmin_m = y.min() - 0.5
    latmax_m = y.max() + 0.5
    lonmin_m = x.min() - 0.5
    lonmax_m = x.max() + 0.5
    #  rango de coordenadas del area a interpolar
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
    z2 = []
    for i, value in enumerate(z):
        if np.isnan(z[i]) is False:
            x2.append(x[i])
            y2.append(y[i])
            z2.append(z[i])
    print(len(x), len(y), len(z))
    # Grilla de interpolacion
    resolucion = 300
    lats = np.linspace(latmin, latmax, resolucion)
    lons = np.linspace(lonmin, lonmax, resolucion)
    y_mapa, x_mapa = np.meshgrid(lats, lons)
    z_mapa = griddata((x, y), z, (x_mapa, y_mapa), method='linear')

    # mask a valores nan
    z_mapa = np.ma.masked_array(z_mapa, mask=np.isnan(z_mapa))

    # plot grilla
    cmap = plt.cm.gnuplot2
    cmap.set_under(color="white", alpha="0")
    # escala para colobar
    z_mapa_esc = z_mapa*1e7  # para arreglar escala
    im = ax.pcolormesh(x_mapa, y_mapa, z_mapa_esc,
                       cmap=cmap,
                       transform=ccrs.PlateCarree(),
                       vmin=z_mapa_esc.min(), vmax=z_mapa_esc.max(),
                       alpha=1, zorder=1)
    cb = plt.colorbar(im)
    cb.set_label('1/year (*1e-7)')

    # detalles finales
    ax.coastlines(resolution='10m')
    return fig
