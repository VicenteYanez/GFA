#! /usr/bin/python3

"""
@author: Vicente Yáñez
@date: 2017
"""

import click

import gfa.scripts.ts_select
import gfa.scripts.ts_build as ts_build
import gfa.scripts.ts_plot as ts_plot
import gfa.scripts.ts_vector as ts_vector
import gfa.scripts.ts_buildmodelall as ts_buildmodelall
import gfa.scripts.ui_buildmodel as ui_buildmodel
import gfa.scripts.gfa_configure as gfa_configure


@click.group()
def cli():
    pass


@cli.command()
def configure():
    """
    Function that creates a default parameters.ini file
    in the current directory
    """
    gfa_configure.main()
    return


@cli.command()
@click.option('--alias', default='default', help='identifier for your query')
@click.option('--lonmin', default=-180., help='Minimun longitude')
@click.option('--lonmax', default=180., help='Maximun longitude')
@click.option('--latmin', default=-90., help='Minimun latitude')
@click.option('--latmax', default=90., help='Minimun latitude')
@click.option('--tmin', default=1., help='Initial time')
@click.option('--tmax', default=3000., help='Final time')
def select(alias, lonmin, lonmax, latmin, latmax, tmin, tmax):
    """
    Extract the time series data giving a certain time range and area
    """
    gfa.scripts.ts_select.main(alias, lonmin, lonmax, latmin, latmax, tmin,
                               tmax)
    return


@cli.command()
@click.option('--alias', default='default', help='identifier for your query')
@click.option('--station', default="", help='Station name')
@click.option('--poly', default=1, help='Grade of the polynomial function')
@click.option('--fourier', default='1',
              help='periods of the Fourier function')
@click.option('--jumps', default='2010.01,2014.24832', help='Step location')
@click.option('--logscale', default='1', help='Scale of the log function')
@click.option('--logstart', default='2010.01',
              help='Start of the log function')
def build(alias, station, poly, jumps, fourier, logscale, logstart):
    """
    Calculates the model for one station
    """
    fourier = [float(s) for s in fourier.split(',')]
    jumps = [float(s) for s in jumps.split(',')]
    logscale = [float(s) for s in logscale.split(',')]
    logstart = [float(s) for s in logstart.split(',')]
    ts_build.main(alias, station, poly, jumps, fourier, logscale, logstart)
    return


@cli.command()
@click.option('--alias', default='', help='identifier for your query')
@click.option('--station', default="", help='Station name')
@click.option('--model', default="", help='-m for model plot')
@click.option('--vector', default="", help='-v for vector plot')
def plot(alias, station, model, vector):
    """
    Plot the time series and optionally include the model and vectors
    """
    ts_plot.main(alias, station, model, vector)
    return


@cli.command()
@click.option('--alias', default='default', help='identifier for your query')
@click.option('--station', default="", help='Station name')
@click.option('--type', default="tangent", help='Vector type, tangent or fit')
@click.option('--period', default=[], help='Time or period of time, it has \
to be a list')
def vector(alias, station, type, period):
    """
    Script that calculates the displacement vector of one station
    in a certain time or period.
    """
    period = [float(s) for s in period.split(',')]
    ts_vector.main(alias, station, type, period)
    return


@cli.command()
def buildmodelall():
    """
    Script that calculates a automatic model for all the station
    """
    ts_buildmodelall.main()
    return


@cli.command()
@click.option('--alias', default='default', help='identifier for your query')
@click.option('--station', default="", help='Station name')
def uimodel(alias, station):
    """
    Runs a manual script that calculates models with the user's parameters.
    """
    ui_buildmodel.main(alias, station)
    return
