
# A very simple Flask Hello World app for you to get started with...

import os
from time import gmtime, strftime
import zipfile

from flask import Flask, render_template, flash, redirect, request, url_for, send_file
import numpy as np

from gfa.log_config import Logger
from gfa.MiddleLayer import MiddleLayer, load_model, load_vectors
from gfa.load_param import Config


filedir = os.path.dirname(os.path.realpath(__file__))

app = Flask(__name__, static_folder='{}/../static'.format(filedir),
            template_folder='{}/../templates'.format(filedir))
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
            latlon, df, sta_param = load_model(model_list_file)

            if os.path.isfile(vector_file):
                # load the vector data if there is a vector file
                df = load_vectors(df, vector_file)
        else:
            df = []
            latlon = []

        return render_template("index.html", stations=np.array(df), latlon=latlon,
                               herenow=strftime("%Y-%m-%d %H:%M:%S", gmtime()))

    except TypeError as e:
        error = 'Exception. Please, check your query. If the error repeats \
contact the admin'
        log = Logger()
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
        log = Logger()
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
        log = Logger()
        log.logger.error(e)
        flash('Sorry, error ploting {}'.format(station))
        return 'error'


@app.route('/vector/', methods=['POST', 'GET'])
def vector():
    """
    Calculate velocity vectors for the station 'station'
    """
    try:
        if request.method == 'POST' and request.form['btn_vector'] == 'Calculate':
            tv1 = request.form['t1']
            tv2 = request.form['t2']
            station = request.form['stationname']

            middle = MiddleLayer(username)
            middle.middle_vector(station, tv1, tv2)
            flash('{} vector sucessfull'.format(station))
    except ValueError as e:
        log = Logger()
        log.logger.error(e)
        flash('Start time should be lower than the end time')
    return redirect(url_for('homepage'))


@app.route('/download')
def download():
    """
    zip the files in the user directory
    """
    output_dir = Config.config['PATH']['output_dir']
    os.chdir(output_dir)
    # check if zipfile exists
    if os.path.isfile('{}/yourcooldata.zip'.format(username)):
        os.remove('{}/yourcooldata.zip'.format(username))

    try:
        zf = zipfile.ZipFile('{}/yourcooldata.zip'.format(username),
                             mode='w')
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
        return send_file('{}{}/yourcooldata.zip'.format(output_dir,username),
                         as_attachment=True)

    except FileNotFoundError as e:
        log = Logger()
        log.logger.error(e)
        flash('error retreaving the zip file')
        return redirect(url_for('homepage'))


# run function
def run():
    app.run()


if __name__ == "__main__":
    app.run()
