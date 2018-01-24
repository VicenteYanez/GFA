
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import numpy as np

from gfa.website_tools import contour2geojson
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
    # ax.coastlines(resolution='10m')
    ax.outline_patch.set_visible(False)
    ax.background_patch.set_visible(False)
    return fig, cbfig


def wk_figure(x, y, z):
    """

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
    zmax = z_mapa.max()
    cmap = LinearSegmentedColormap.from_list('mycmap',
                                             [(0, 'blue'),
                                              (0.8/zmax, 'turquoise'),
                                              (1.2/zmax, 'yellow'),
                                              (1, 'red')], N=256)
    # colorbar scale
    im = ax.pcolormesh(x_mapa, y_mapa, z_mapa, vmin=z_mapa.min(), vmax=zmax,
                       cmap=cmap, transform=ccrs.PlateCarree(), alpha=1, zorder=1)

    cbfig = plt.figure()
    cbax = plt.axes(projection=ccrs.Mercator())
    cbar = cbfig.colorbar(im)
    # cbar.set_ticks(self.config.colorbar_ticks)
    cbax.set_visible(False)

    # final details
    # ax.coastlines(resolution='10m')
    ax.outline_patch.set_visible(False)
    ax.background_patch.set_visible(False)
    return fig, cbfig


def s_figure(x, y, z):
    """

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
    zmax = z_mapa.max()
    cmap = plt.cm.Reds
    # colorbar scale
    im = ax.pcolormesh(x_mapa, y_mapa, z_mapa, vmin=z_mapa.min(), vmax=zmax,
                       cmap=cmap, transform=ccrs.PlateCarree(), alpha=1,
                       zorder=1)

    cbfig = plt.figure()
    cbax = plt.axes(projection=ccrs.Mercator())
    cbar = cbfig.colorbar(im)
    cbar.set_label('|S|')
    # cbar.set_ticks(self.config.colorbar_ticks)
    cbax.set_visible(False)

    # final details
    # ax.coastlines(resolution='10m')
    ax.outline_patch.set_visible(False)
    ax.background_patch.set_visible(False)
    return fig, cbfig


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
                                             alfa=alfa, dmin=100000,
                                             method='numpy')

    # tensor calculations
    S, W = field.velocitytensor_2d(gradiente[0], gradiente[1],
                                   gradiente[2], gradiente[3])
    wz = field.vertical_vorticity2d(W)
    Wk = field.cinematic_vorticity(S, W)
    S2 = field.frobenius_norm(S)

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
    cbfig.savefig("{}/wz_colorbar.png".format(dir_path), dpi=150,
                  bbox_inches='tight', pad_inches=0, transparent=True)
    fname = "{}/wz_field.png".format(dir_path)
    fig1.savefig(fname, dpi=500, edgecolor='k', bbox_inches='tight',
                 transparent=True, pad_inches=0.)
    # cinematic vorticity figure
    fig2, cbfig_wk = wk_figure(gridx, gridy, Wk)
    cbfig_wk.savefig("{}/wk_colorbar.png".format(dir_path), dpi=150,
                     bbox_inches='tight', pad_inches=0, transparent=True)
    fname2 = "{}/wk_field.png".format(dir_path)
    fig2.savefig(fname2, dpi=500, edgecolor='k', bbox_inches='tight',
                 transparent=True, pad_inches=0.)
    # cinematic vorticity figure
    fig3, cbfig_s = s_figure(gridx, gridy, S2)
    cbfig_s.savefig("{}/s2_colorbar.png".format(dir_path), dpi=150,
                    bbox_inches='tight', pad_inches=0, transparent=True)
    fname3 = "{}/s2_field.png".format(dir_path)
    fig3.savefig(fname3, dpi=500, edgecolor='k', bbox_inches='tight',
                 transparent=True, pad_inches=0.)

    print('end function')
    return


def vectors_map(dir_path, x, y, ve, vn, lat_range, lon_range):
    # plot histogram of values
    f, (ax1, ax2) = plt.subplots(2, figsize=(12, 8))
    ax1.hist(ve, bins=30)
    ax2.hist(vn, bins=30)
    f.savefig("{}/vector_hist.png".format(dir_path))
    plt.close()

    fig = plt.figure()
    ax = plt.axes(projection=ccrs.Mercator())
    ax.outline_patch.set_visible(False)
    ax.background_patch.set_visible(False)

    # zone to interpol
    latmin = lat_range[0]
    latmax = lat_range[1]
    lonmin = lon_range[0]
    lonmax = lon_range[1]

    extend = [lonmin, lonmax, latmin, latmax]
    ax.set_extent(extend)

    vector = plt.quiver(x, y, ve, vn, transform=ccrs.PlateCarree())
    plt.quiverkey(vector, 0., 1.015, 0.03, u'30 mm/año', labelpos="E")
    fname = "{}/vectorsmap.png".format(dir_path)
    fig.savefig(fname, dpi=500, bbox_inches='tight', transparent=True,
                pad_inches=0.)
    plt.close()

    # ax.background_img(name="ne_ultra", resolution='uhigh', extent=extend)
    scalefig = plt.figure()
    scaleax = plt.axes(projection=ccrs.Mercator())
    # cbar.set_ticks(self.config.colorbar_ticks)
    scaleax.set_visible(False)
    plt.quiverkey(vector, 0., 1.015, 0.03, u'30 mm/año', labelpos="E")
    scalefig.savefig("{}/vectorscale.png".format(dir_path), dpi=150,
                     bbox_inches='tight', pad_inches=0, transparent=True)
    plt.close()

    return
