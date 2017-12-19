
from datetime import datetime
import time
import os
import shutil
import json
import re

import numpy as np
from flask import flash
import pandas as pd

import gfa.log_config as log
from gfa.load_param import Config
from gfa.field_analysis.cardozo_vorticity import cardozo_vorticity
from gfa.gnss_analysis.VectorData import VectorData


"""
This is a ugly file
The problem is that GFA loads the parameter.ini file that it is
inside the CURRENT WORKING DIRECTORY so before import gfa.load_param
i need to change it.
I hope i'll be fixing this in the next version of GFA.
"""
output_dir = Config.config['PATH']['output_dir']
os.chdir(output_dir)


def toYearFraction(date):
    def sinceEpoch(date):  # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction


class MiddleLayer():
    """
    Link class between flask app and GFA
    """
    def __init__(self, user):
        from gfa.load_param import Config
        self.user = user
        self.lista_gps = Config.config['PATH']['ListaGPS']
        self.series_dir = Config.config['PATH']['GPSdata']
        self.output_dir = Config.config['PATH']['output_dir']
        self.user_dir = "{}{}".format(self.output_dir, self.user)

    def middle_select(self, lonmin, lonmax, latmin, latmax, tmin, tmax):
        if os.path.isdir(self.user_dir):
            # remove previus selection
            shutil.rmtree(self.user_dir)
        ti = datetime.strptime(tmin, "%Y-%m-%d")
        tf = datetime.strptime(tmax, "%Y-%m-%d")

        # time format transform
        ti = float(toYearFraction(ti))
        tf = float(toYearFraction(tf))

        try:
            # selection
            from gfa.scripts import ts_select
            ts_select.main(self.user, lonmin, lonmax, latmin, latmax, ti,
                           tf)
            # save select param
            selectfile = '{}/select_param.json'.format(self.user_dir)
            paramdict = {
                'Select longitude min': lonmin,
                'Select longitude max': lonmax,
                'Select latitude min': latmin,
                'Select latitude max': latmax,
                'Select time start': ti,
                'Select time end': tf
            }

            if os.path.isfile(selectfile):
                os.remove(selectfile)

            with open(selectfile, 'xt') as f:
                json.dump(paramdict, f)
        except (ImportError, TypeError, ValueError) as e:
            error = 'Error importing or working with GFA. {}'.format(e)
            flash('Exception. Please, check your query. If the error repeats \
contact the admin')
            log.logger.error(error)
            return []

        return

    def middle_edit(self, station, poly, jumps, fourier, logstart, logscale):
        try:
            # selection
            from gfa.scripts import ts_build
            ts_build.main(self.user, station, poly, jumps, fourier, logscale,
                          logstart)
        except (ImportError) as e:
            error = 'Error importing or working with GFA. {}'.format(e)
            flash('Exception. Please, check your query. If the error repeats \
contact the admin')
            log.logger.error(error)
        return

    def middle_plot(self, station, vector='-'):
        """
        vector parameter is optional, if is given should and
        is equal to '-v' string, the vectors are plot
        """
        # add verification if the file exist
        output_plot = "{}/{}.png".format(self.user_dir, station)
        if os.path.isfile(output_plot):
            os.remove(output_plot)

        from gfa.scripts import ts_plot
        if vector == '-v':
            ts_plot.main(self.user, station, '-m', '-v')
        else:
            ts_plot.main(self.user, station, '-m', '-')

        return output_plot

    def middle_vector(self, station, ti, tf, vector_file):
        """

        """
        # ti and tf are a string with parenthesis, first of all
        # we need to remove ir and turning they in a
        # list of floats
        ti = [float(toYearFraction(datetime.strptime(s, "%Y-%m-%d")))
              for s in ti.split(',')]
        tf = [float(toYearFraction(datetime.strptime(s, "%Y-%m-%d")))
              for s in tf.split(',')]
        if os.path.isfile(vector_file):
            vdata = VectorData(vector_file)
            tif, remove_times = vdata.check_time(station, ti, tf)
        else:
            tif = np.array([ti, tf]).T
            remove_times = []

        # calculate and save the vectors
        from gfa.scripts import ts_vector
        for times in tif:
            # check value of time
            # parameter vtype is for ts_vector() function
            if times[0] > times[1]:
                raise ValueError
            elif times[0] == times[1]:
                vtype = 'tangent'
            elif times[0] < times[1]:
                vtype = 'fit'
            ts_vector.main(self.user, station, vtype, [times[0], times[1]])

        # remove the user erased vectors
        for times in remove_times:
            vdata.remove_vector(station, remove_times)
        return

    def middle_field(self, model_list_file, vector_file, lon_range, lat_range,
                     time_range, grid, alfa):
        # read date
        ti = datetime.strptime(time_range[0], "%Y-%m-%d")
        tf = datetime.strptime(time_range[1], "%Y-%m-%d")

        # load df with the position of the station
        latlon, df = load_model(model_list_file)

        # load vectors from file
        if os.path.isfile(vector_file):
            vdata = VectorData(vector_file)
            x, y, ve, vn = vdata.select(df, lon_range, lat_range, [ti, tf])
            pngfile = "{}/field.png".format(self.user_dir)
            if os.path.isfile(pngfile):
                os.remove(pngfile)
            cardozo_vorticity(self.user_dir, x, y, ve, vn,
                              lat_range, lon_range, grid, alfa)

            # save select param
            fieldfile = '{}/select_param.json'.format(self.user_dir)
            paramdict = {
                'Field longitude min': lon_range[0],
                'Field longitude max': lon_range[1],
                'Field latitude min': lat_range[0],
                'Field latitude max': lat_range[1],
                'Field time start': time_range[0],
                'Field time end': time_range[1]
            }

            with open(fieldfile, 'wt') as f:
                json.dump(paramdict, f)

        else:
            flash('no vector file, please calculate some vectors first!')

        return

    def middle_vectortable(self, df, vector_file):
        """
        """
        vdata = VectorData(vector_file)
        vector_table = vdata.get_table(df)

        return vector_table


def load_model(model_list_file):
    """
    Load the model list file
    """
    # loads name of stations and json with the parameters
    data = np.loadtxt(model_list_file, delimiter="    ",
                      skiprows=1, dtype=bytes).astype(str).T
    # json parameters to dictionary
    parameters = [json.loads(parameter) for parameter in data[3]]
    poly = [p["polinomio"] for p in parameters]
    jumps = [j["saltos"] for j in parameters]
    fourier = [f["Periodos Fourier"] for f in parameters]
    log_i = [l["Inicio log"] for l in parameters]
    log_sc = [s["Escala curva log"] for s in parameters]

    # list with the station name and its parameters
    sta_param = [data[0], data[1], data[2], poly, jumps, fourier,
                 log_i, log_sc]
    sta_param = [list(x) for x in zip(*sta_param)]  # transpose data
    latlon = np.array([data[1], data[2]]).T.astype(float)

    df = pd.DataFrame(sta_param, columns=['station', 'longitude',
                                          'latitude', 'polynomial',
                                          'jumps', 'fourier',
                                          'log start', 'log scale'])
    df['longitude'] = df['longitude'].apply(
        lambda x: '{0:.3f}'.format(float(x)))
    df['latitude'] = df['latitude'].apply(lambda x: '{0:.3f}'.format(float(x)))
    return latlon, df
