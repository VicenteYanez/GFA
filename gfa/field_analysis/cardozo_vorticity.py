
import copy

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import numpy as np

from gfa.field_analysis import contour2geojson
from gfa.field_analysis import field


def greater_absolute(value1, value2):
    if abs(value1) >= abs(value2):
        return np.array([-value1, value1])
    else:
        return np.array([value2, -value2])


def data_stats(lon, lat, ve, vn):

    return


def wz_figure(x, y, z):
    """
    Figure of example_vorticity.py
    """
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.Mercator())

    #  zone to interpol
    latmin = y.min()
    latmax = y.max()
    lonmin = x.min()
    lonmax = x.max()

    extend = [lonmin, lonmax, latmin, latmax]
    ax.set_extent(extend)

    # ax.background_img(name="ne_ultra", resolution='uhigh', extent=extend)

    # interpolation grid
    resolucion = 1000
    lats = np.linspace(latmin, latmax, resolucion)
    lons = np.linspace(lonmin, lonmax, resolucion)
    y_mapa, x_mapa = np.meshgrid(lats, lons)
    z_mapa = griddata((x, y), z, (x_mapa, y_mapa), method='linear')

    # mask nan values
    z_mapa = np.ma.masked_array(z_mapa, mask=np.isnan(z_mapa))

    # plot grid
    cmap = plt.cm.coolwarm

    # colorbar scale
    vesc = greater_absolute(z_mapa.min(), z_mapa.max())
    im = ax.pcolormesh(x_mapa, y_mapa, z_mapa,
                       cmap=cmap, transform=ccrs.PlateCarree(),
                       vmin=vesc.min(), vmax=vesc.max(), alpha=1, zorder=1)

    cbfig = plt.figure()
    cbax = plt.axes(projection=ccrs.Mercator())
    cbar = cbfig.colorbar(im)
    cbar.set_label('vertical vorticity [degrees/year]')
    # cbar.set_ticks(self.config.colorbar_ticks)
    cbax.set_visible(False)

    # final details
    ax.coastlines(resolution='10m')
    return fig, cbfig


def wk_figure(x, y, z):
    """
    Figure of example_vorticity.py
    """
    fig = plt.figure(figsize=(9, 19))
    ax = plt.axes(projection=ccrs.PlateCarree())

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

    paralelos = np.arange(int(latmin_m), int(latmax_m), 2)
    meridianos = np.arange(int(lonmin_m), int(lonmax_m), 2)

    extend = [lonmin_m, lonmax_m, latmin_m, latmax_m]
    ax.set_extent(extend)
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

    ax.background_img(name="ne_ultra", resolution='uhigh', extent=extend)

    # interpolation grid
    resolucion = 400
    lats = np.linspace(latmin, latmax, resolucion)
    lons = np.linspace(lonmin, lonmax, resolucion)
    y_mapa, x_mapa = np.meshgrid(lats, lons)
    z_mapa = griddata((x, y), z, (x_mapa, y_mapa), method='linear')

    # mask nan values
    z_mapa = np.ma.masked_array(z_mapa, mask=np.isnan(z_mapa))
    zmin = z_mapa.min()
    zmax = z_mapa.max()

    # plot grid
    cmap = LinearSegmentedColormap.from_list('mycmap',
                                             [(0, 'blue'),
                                              (0.8/zmax, 'turquoise'),
                                              (1.2/zmax, 'yellow'),
                                              (1, 'red')], N=256)
    # colorbar scale
    im = ax.pcolormesh(x_mapa, y_mapa, z_mapa,
                       cmap=cmap,
                       transform=ccrs.PlateCarree(),
                       vmin=zmin, vmax=zmax,
                       alpha=0.5, zorder=1)
    plt.colorbar(im)

    # final details
    ax.coastlines(resolution='10m')
    return fig


def residual_figure(x, y, z):
    return


def fun(a, M):
    return M*a


def vorticity_contour(dir_path, x, y, z):
    """
    Function for create a contourn and transform it
    to a geojson file
    """

    figure = plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    # mask nan values
    z = np.ma.masked_array(z, mask=np.isnan(z))
    # geojsoncontour code
    n_contours = 15
    levels = np.linspace(start=z.min(), stop=z.max(), num=n_contours)

    Z = np.reshape(z, (len(x), len(x[0])))
    contour = ax.contour(x, y, Z, levels=levels, cmap=plt.cm.coolwarm)

    # save file
    fname = "{}/out.geojson".format(dir_path)

    # Convert matplotlib contour to geojson
    contour2geojson.contour_to_geojson(
        contour=contour,
        geojson_filepath=fname,
        contour_levels=levels,
        min_angle_deg=10.0,
        ndigits=3,
        unit='degrees/year'
    )
    return


def cardozo_vorticity(dir_path, lon_gps, lat_gps, ve_gps, vn_gps,  lat_range,
                      lon_range, grid, alfa=50000):
    # make grid
    gridx2d, gridy2d = np.meshgrid(np.linspace(lon_range[0], lon_range[1], 20),
                               np.linspace(lat_range[0], lat_range[1], 60))
    gridx = gridx2d.ravel()
    gridy = gridy2d.ravel()

    # conversion mm/año a m/año
    ve_gps = ve_gps/1000
    vn_gps = vn_gps/1000

    # cardozo&allmendinger strain
    b = np.reshape([ve_gps, vn_gps], (2*len(ve_gps), 1), order='F')
    gradiente, r = field.distance_weigthed2d(b, lon_gps, lat_gps, gridx, gridy,
                                             alfa=alfa, dmin=100000)

    # tensor calculations
    S, W = field.velocitytensor_2d(gradiente[0], gradiente[1],
                                   gradiente[2], gradiente[3])
    wz = field.vertical_vorticity2d(W)
    Wk = field.cinematic_vorticity(S, W)

    # save data
    wz = wz*360  # to convert to degrees/year
    savefilewz = '{}/wz.txt'.format(dir_path)
    savefilewk = '{}/wk.txt'.format(dir_path)
    np.savetxt(savefilewz, np.array([gridx, gridy, wz]).T,
               header='x y z(degrees/year)')
    np.savetxt(savefilewk, np.array([gridx, gridy, Wk]).T, header='x y wk')

    # plot histogram of values
    f, (ax1, ax2) = plt.subplots(2, figsize=(12, 8))
    ax1.hist(np.ma.masked_array(wz, mask=np.isnan(wz)), bins=30)
    ax2.hist(np.ma.masked_array(Wk, mask=np.isnan(Wk)), bins=30)
    f.savefig("{}/hist.png".format(dir_path))

    # Graficar
    # vertical vorticity figure
    vorticity_contour(dir_path, gridx2d, gridy2d, wz)
    fig1, cbfig = wz_figure(gridx, gridy, wz)
    cbfig.savefig("{}/colorbar.png".format(dir_path), dpi=150,
                  bbox_inches='tight', pad_inches=0, transparent=True)
    fname = "{}/wz_field.png".format(dir_path)
    fig1.savefig(fname, dpi=500, edgecolor='k', bbox_inches='tight',
                 transparent=True, pad_inches=0.)
    # cinematic vorticity figure
    fig2 = wk_figure(gridx, gridy, Wk)
    fname2 = "{}/wk_field.png".format(dir_path, alfa)
    fig2.savefig(fname2, dpi=300, edgecolor='k', orientation='portrait',
                 bbox_inches=None, pad_inches=0.2)

    print('end function')
    return
