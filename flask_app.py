
# A very simple Flask Hello World app for you to get started with...

import json
import os
from time import gmtime, strftime
import zipfile

from flask import Flask, render_template, flash, redirect, request, url_for, send_file
import numpy as np

import gfa.log_config as log
from MiddleLayer import MiddleLayer
from gfa.load_param import Config


app = Flask(__name__)
app.secret_key = 'super secret key'
app.config['SESSION_TYPE'] = 'filesystem'

userdir = './gfa'
username = 'gfa'


@app.errorhandler(404)
def page_not_found(e):
    return render_template("404.html")


@app.errorhandler(405)
def method_not_found(e):
    return render_template("405.html")


@app.route('/', methods=['GET', 'POST'])
def homepage():
    model_list_file = "{}/modelo_lista.txt".format(userdir)
    vector_file = "{}/vectors.txt".format(userdir)
    try:
        if request.method == 'POST' and request.form['btn_select'] == 'Select':
            latmin = float(request.form['latmin'])
            latmax = float(request.form['latmax'])
            lonmin = float(request.form['lonmin'])
            lonmax = float(request.form['lonmax'])
            ti = request.form['t_inicial']
            tf = request.form['t_final']

            middle = MiddleLayer(username)
            middle.middle_select(lonmin, lonmax, latmin, latmax, ti, tf)
            flash('finished the data loading')

        if os.path.isfile(model_list_file):
            # if there is data, load it
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
            alldata = [data[0], data[1], data[2], poly, jumps, fourier,
                       log_i, log_sc]
            alldata = [list(x) for x in zip(*alldata)]  # transpose data
            latlon = np.array([data[1], data[2]]).T.astype(float)

            if os.path.isfile(vector_file):
                # load the vector data if there is a vector file
                vectordata = np.loadtxt(vector_file, delimiter="    ",
                                        usecols=[0, 2, 3, 4, 8, 9],
                                        skiprows=1, dtype=bytes).astype(str)
                #
                # join vector data with the rest of the data
                vectordata = [v[0], v[1]]

        else:
            alldata = []

        return render_template("index.html", stations=alldata, latlon=latlon,
                               herenow=strftime("%Y-%m-%d %H:%M:%S", gmtime()))

    except TypeError as e:
        error = 'Exception. Please, check your query. If the error repeats \
contact the admin'
        log.logger.error(e)
        flash(error)
    return render_template("index.html", stations=[])


@app.route('/about/')
def aboutpage():
    return render_template('about.html')


@app.route('/edit/<station>', methods=['POST', 'GET'])
def edit(station):
    """
    Code for handle edit request
    """
    # transform input to a list with floats
    new_poly = int(request.form['poly'])
    new_jump = list(np.array([x for x in request.form['jump'].strip('[]').split(',') if x != '']).astype(float))
    new_fourier = list(np.array([x for x in request.form['fourier'].strip('[]').split(',') if x != '']).astype(float))
    new_logstart = list(np.array([x for x in request.form['logstart'].strip('[]').split(',') if x != '']).astype(float))
    new_logscale = list(np.array([x for x in request.form['logscale'].strip('[]').split(',') if x != '']).astype(float))

    middle = MiddleLayer(username)
    middle.middle_edit(station, new_poly, new_jump, new_fourier, new_logstart,
                       new_logscale)

    flash("{} upgrated".format(station))

    return redirect(url_for('homepage'))


@app.route('/plot/<station>')
def plots(station):
    """
    Plot request handlign
    """
    try:
        middle = MiddleLayer(username)
        filename = middle.middle_plot(station)
        # cache_timeout should stop that t
        return send_file(filename)
    except (ValueError, IOError) as e:
        log.logger.error(e)
        flash('Sorry, error ploting {}'.format(station))
        return 'error'


@app.route('/plot-v/<station>')
def plots_v(station):
    """
    Plot request handlign
    """
    try:
        middle = MiddleLayer(username)
        filename = middle.middle_plot(station, '-v')
        # cache_timeout should stop that t
        return send_file(filename)
    except (ValueError, IOError) as e:
        log.logger.error(e)
        flash('Sorry, error ploting {}'.format(station))
        return 'error'


@app.route('/vector/<station>')
def vector(station):
    """
    Calculate velocity vectors for the station 'station'
    """
    return


@app.route('/download/')
def download():
    """
    zip the files in the user directory
    """
    output_dir = Config.config['PATH']['output_dir']
    os.chdir(output_dir)
    zf = zipfile.ZipFile('{}/yourcooldata.zip'.format(username),
                         mode='w')
    try:
        zf.write('{}/series_lista.txt'.format(username), 'list_series.txt')
        zf.write('{}/modelo_lista.txt'.format(username), 'list_model.txt')

        for seriesfile in os.listdir('{}/series'.format(username)):
            zf.write('{}/series/{}'.format(username, seriesfile),
                     'series/{}'.format(seriesfile))
        for modelfile in os.listdir('{}/modelo'.format(username)):
            zf.write('{}/modelo/{}'.format(username, modelfile),
                     'model/{}'.format(modelfile))
        for pngfile in [png for png in os.listdir(username)
                        if png.endswith('.png')]:
            zf.write('{}/{}'.format(username, pngfile), pngfile)
        vectorfile = '{}/vectors.txt'.format(username)
        if os.path.isfile(vectorfile):
            zf.write(vectorfile, vectorfile)
        zf.close()
        return send_from_directory('{}{}'.format(output_dir, username),
                                   'yourcooldata{}.zip'.format(strftime("%Y-%m-%d %H:%M:%S", gmtime())),
                                   as_attachment=True)

    except FileNotFoundError as e:
        print(e)
        flash('error retreaving the zip file')
        return redirect(url_for('homepage'))


if __name__ == "__main__":
    app.run()
