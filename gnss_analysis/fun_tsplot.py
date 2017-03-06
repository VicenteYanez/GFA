#! /usr/bin/python3

"""
@author: Vicente Yáñez

Functions for plot time series GNSS
"""

import pdb

import numpy as np
import matplotlib.pyplot as plt


def plot(estac, tiempo, desp):
    """
    Grafica o guarda las figuras construidas.
    Input
    estac   : lista de string con nombre estaciones
    tiempo  : array 1d
    desp    : array de 2d de 3 x n. Donde cada fila es el desp en un eje.
              [0] Este, [1] Norte [3] Vertical
    modelo  : Idem que desp
    """
    # plot puntos y error bar
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize=(10, 7))

    # formato ejes
    plt.xlabel('years', fontsize=8)
    plt.ylabel('mm', fontsize=8, position=(0.5, 1.8))

    # change font size x and y axis
    plt.xticks(fontsize=8)
    ax1.tick_params(labelsize=8)
    ax2.tick_params(labelsize=8)
    ax3.tick_params(labelsize=8)

    # poner marcas solo en parte inferior y laterales
    ax1.get_xaxis().tick_bottom()
    ax2.get_xaxis().tick_bottom()
    ax3.get_xaxis().tick_bottom()

    # zero line
    ax1.axhline(0, color='#9D9D9D', linestyle='--', linewidth=0.5)
    ax2.axhline(0, color='#9D9D9D', linestyle='--', linewidth=0.5)
    ax3.axhline(0, color='#9D9D9D', linestyle='--', linewidth=0.5)

    # evitar offser
    ax1.ticklabel_format(useOffset=False)
    ax2.ticklabel_format(useOffset=False)
    ax3.ticklabel_format(useOffset=False)

    # title
    ax1.set_title("{} Este".format(estac), fontsize=9)
    ax2.set_title("{} Norte".format(estac), fontsize=9)
    ax3.set_title("{} Vertical".format(estac), fontsize=9)

    # plot data
    ax1.plot(tiempo, desp[0], 'o', color='#99A3F2', markeredgewidth=0.0)
    ax2.plot(tiempo, desp[1], 'o', color='#99A3F2', markeredgewidth=0.0)
    ax3.plot(tiempo, desp[2], 'o', color='#99A3F2', markeredgewidth=0.0)

    return f, (ax1, ax2, ax3)


def add_modelo(figure, axes, tiempo, modelo):
    axes[0].plot(tiempo, modelo[0], '-g', label='Este')
    axes[1].plot(tiempo, modelo[1], '-g', label='Norte')
    axes[2].plot(tiempo, modelo[2], '-g', label='Vertical')

    return figure, axes


def add_velocity(figure, axes, tiempo, rectas):
    """
    Add velocity trend
    """
    tpart = (tiempo[-1] - tiempo[0])/5
    # muestreo
    for recta in rectas:
        if recta[2][0] == recta[2][1]:
            t_tangent = recta[2][0]
            tr = [t for t in tiempo if t > (t_tangent-tpart) and
                  t < (t_tangent+tpart)]
        else:
            tr = [t for t in tiempo if t > recta[2][0] and t < recta[2][1]]

        d1 = recta[0]*tr[0] + recta[1]
        d2 = recta[0]*tr[-1] + recta[1]
        axes[0].plot((tr[0], tr[-1]), (d1[0], d2[0]),
                     '-r', linewidth=0.8, label='Velocidad Este')
        axes[1].plot((tr[0], tr[-1]), (d1[1], d2[1]),
                     '-r', linewidth=0.8, label='Velocidad Norte')
        axes[2].plot((tr[0], tr[-1]), (d1[2], d2[2]),
                     '-r', linewidth=0.8, label='Velocidad Vertical')
    return figure, axes


def add_sismo(figure, axes, t_sismo):
    # agregar evento sismico
    for sismo in t_sismo:
        axes[0].plot((sismo, sismo), (desp[0].max(), desp[0].min()), '-k')
        axes[1].plot((sismo, sismo), (desp[1].max(), desp[1].min()), '-k')
        axes[2].plot((sismo, sismo), (desp[2].max(), desp[2].min()), '-k')
    return figure, axes


def save_figure(figure, savefile='./',):
    # Guardar figura
    plt.savefig(savefile, dpi=300, bbox_inches='tight')
    plt.close()


def mapa_estaciones(lat, lon, pathdir='./'):
    """
    Grafica la ubicacion de estaciones como puntos en un mapa
    """
    # rango de coordenadas del mapa
    latmin = np.array(lat).min()-.5
    latmax = np.array(lat).max()+.5
    lonmin = np.array(lon).min()-.5
    lonmax = np.array(lon).max()+.5

    paralelos = np.arange(latmin, latmax, 10)
    meridianos = np.arange(lonmin, lonmax, 5)

    proj = ccrs.PlateCarree()

    ax = plt.axes(projection=proj)
    ax.set_extent([lonmin, lonmax, latmin, latmax])
    ax.set_xticks(meridianos, crs=proj)
    ax.set_yticks(paralelos, crs=proj)
    ax.gridlines(crs=proj)

    map.scatter(lon, lat, c="b", transform=proj)
    ax.coastlines(resolution='10m')

    plt.savefig(self.workspace+'mapa_estaciones.ps',
                dpi=300, bbox_inches='tight')
    plt.close()

    return


def add_jlabels():
    """
    Function that add a label where the jump is located
    """
    return


def mapa_vectores(lat, lon, vel_e, vel_n, pathdir='./', filename='vector_map'):
    """
    Grafica un mapa con la velocidad de estaciones dadas
    """
    # rango de coordenadas del mapa
    latmin = np.array(lat).min()-.5
    latmax = np.array(lat).max()+.5
    lonmin = np.array(lon).min()-.5
    lonmax = np.array(lon).max()-.5

    paralelos = np.arange(-56, -17, 5)
    meridianos = np.arange(-80, -63, 2)

    proj = ccrs.PlateCarree()

    ax = plt.axes(projection=proj)
    ax.set_extent([lonmin, lonmax, latmin, latmax])
    ax.set_xticks(meridianos, crs=proj)
    ax.set_yticks(paralelos, crs=proj)
    ax.gridlines(crs=proj)

    vector = plt.quiver(x, y, vel_e, vel_n, width=0.004, color='k')
    plt.quiverkey(vector, 0.8, 1.02, 25, u'25 mm/año',
                  fontproperties={'size': 10}, labelpos="E",
                  transform=proj)

    plt.savefig('{}{}.jpg'.format(pathdir, filename), dpi=400,
                bbox_inches='tight', pad_inches=0.2)

    plt.close()
    return
