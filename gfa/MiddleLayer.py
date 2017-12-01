
from datetime import datetime
import time
import os
import shutil
import json

import numpy as np
from flask import flash
import pandas as pd

import gfa.log_config as log
from gfa.load_param import Config

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
        self.user_dir = "{}/{}".format(self.output_dir, self.user)

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

    def middle_vector(self, station, ti, tf):
        """

        """

        ti = datetime.strptime(ti, "%Y-%m-%d")
        tf = datetime.strptime(tf, "%Y-%m-%d")

        # time format transform
        ti = float(toYearFraction(ti))
        tf = float(toYearFraction(tf))

        # check value of time
        # parameter vtype is for ts_vector() function
        if ti > tf:
            print(ti, tf)
            raise ValueError
        elif ti == tf:
            vtype = 'tangent'
        elif ti != tf:
            vtype = 'fit'

        from gfa.scripts import ts_vector
        ts_vector.main(self.user, station, vtype, [ti, tf])

        return


def load_vectors(df, vector_file):
    """
    Loads the vector file and join it with the rest of the statio data
    """
    vectordata = np.loadtxt(vector_file, delimiter="    ",
                            usecols=[0, 2, 3, 4, 8, 9],
                            skiprows=1, dtype=bytes).astype(str)
    vector_df = pd.DataFrame(vectordata, columns=['station',
                                                  've', 'vn', 'vz',
                                                  't1', 't2'])
    df = pd.merge(df, vector_df, how='left', on='station')

    return df


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
    return latlon, df, sta_param
