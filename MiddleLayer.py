
from datetime import date, datetime, timedelta
import time
import os
import shutil

import numpy as np
from flask import flash

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

    def middle_vector(self):
        return
