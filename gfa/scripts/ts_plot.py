#! /usr/bin/python3

"""
@author: Vicente Yáñez

Generates a time series figure with the gnss and model data
Parameters
[1] identifier code
[2] station to plot
[3] If the third parameter is -m, besides the time series, the model
is plot.
[4] If the forth parameter is -v, the velocity vector is plot
"""

from gfa.gnss_analysis.loadGPS import load1stations
from gfa.gnss_analysis.model_accions import load_vector
import gfa.gnss_analysis.fun_tsplot as pltfun
from gfa.log_config import Logger
from gfa.load_param import Config


def main(codigo, station, addmodel, addvector):
    # select directory
    directory = '{}{}/'.format(Config.config['PATH']['output_dir'], codigo)
    savefile = '{}{}.png'.format(directory, station)

    # load data and plot function
    tsdata, modeldata = load1stations(station, directory, True, True)
    x0 = tsdata[1][1][0]
    y0 = tsdata[1][2][0]
    z0 = tsdata[1][3][0]
    tsdata[1][1] = tsdata[1][1] - x0
    tsdata[1][2] = tsdata[1][2] - y0
    tsdata[1][3] = tsdata[1][3] - z0
    f, axes = pltfun.plot(station, tsdata[1][0], tsdata[1][1:4])

    # optional plots: add model
    if addmodel == '-m':
        modeldata[1][1] = modeldata[1][1] - x0
        modeldata[1][2] = modeldata[1][2] - y0
        modeldata[1][3] = modeldata[1][3] - z0
        f, axes = pltfun.add_modelo(f, axes, modeldata[1][0],
                                    modeldata[1][1:4])
        print('Model added')
    else:
        print('Model not added')

    # optional plots: add vector
    if addvector == '-v':
        # file vector
        filevector = '{}vectors.txt'.format(directory)
        try:
            # load vectors
            vectors = load_vector(filevector, station)
            # plot
            f, axes = pltfun.add_velocity(f, axes, tsdata[1][0], vectors)
            print('Vector added')
        except FileNotFoundError as err:
            print("Error: Vector file not found")
            log = Logger()
            log.logger.error(err)

    else:
        print('Vector not added')

    # save plot
    pltfun.save_figure(f, savefile)
