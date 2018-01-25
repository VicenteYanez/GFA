
# A very simple Flask Hello World app for you to get started with...
from functools import wraps
import os
from time import gmtime, strftime
import zipfile
import copy
import json
import gc
import sys
import linecache
import shutil

from flask import Flask, render_template, flash, redirect, request, url_for
from flask import session, send_file
import numpy as np

from gfa.log_config import Logger
from gfa.MiddleLayer import MiddleLayer, load_model
from gfa.load_param import Config
from gfa.gnss_analysis.fun_vector import TimeSeriesError


def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)

    msg = 'EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno,
                                                       line.strip(), exc_obj)
    return msg


filedir = os.path.dirname(os.path.realpath(__file__))
# list with the time series databases available
gnss_databases = [s for s in Config.config['WEB']['databases'].split(',')]

app = Flask(__name__, static_folder='{}/../static'.format(filedir),
            template_folder='{}/../templates'.format(filedir))
app.secret_key = 'super secret key'
app.config['SESSION_TYPE'] = 'filesystem'


# if gfa is loaded into a website connect to a db
websession = Config.config['WEB']['websession']
if websession == 'True':
    from flask_sqlalchemy import SQLAlchemy

    # SQLAlchemy config
    # SQLAlchemy declarative method
    SQLALCHEMY_DATABASE_URI = "postgresql://pi:okada1985@localhost/andes3db"
    app.config["SQLALCHEMY_DATABASE_URI"] = SQLALCHEMY_DATABASE_URI
    app.config["SQLAlCHEMY_POOL_RECYCLE"] = 299
    app.config["SQLALCHEM_TRACK_MODIFICATIONS"] = False
    app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

    db = SQLAlchemy(app)

    class User(db.Model):
        """
        Model for the db connection with SQLAlchemy
        """
        __tablename__ = "users"

        ide = db.Column(db.Integer, primary_key=True)
        user = db.Column(db.String(50))
        password = db.Column(db.String(50))

        def __init__(self, user, password):
            self.user = user
            self.password = password

        def __repr__(self):
            return "<User(user={}, password={})>".format(self.user,
                                                         self.password)


def login_required(f):
    @wraps(f)
    def wrap(*args, **kwars):
        if 'logged_in' in session:
            return f(*args, **kwars)
        else:
            return redirect(url_for('login_page'))
    return wrap


@app.route('/login/', methods=['GET', 'POST'])
def login_page():
    error = None
    try:
        if websession == 'False':
            session['logged_in'] = True
            session['username'] = 'gfa'
            return redirect(url_for('homepage'))
        if request.method == 'POST':
            attemped_password = request.form['password']
            attemped_username = request.form['username']
            pwd = User.query.filter_by(user=attemped_username).first().password

            if attemped_password == pwd:
                session['logged_in'] = True
                session['username'] = request.form['username']
                gc.collect()
                return redirect(url_for('homepage'))
            else:
                error = "Invalid credentials. Try again"
                flash(error)
                return render_template("login.html")
        return render_template("login.html")

    except Exception as e:
        message = "A error has been raised, please contact the webmaster"
        log = Logger()
        log.logger.error(PrintException())
        flash(message)
        return render_template("login.html")


@app.errorhandler(404)
def page_not_found(e):
    return render_template("404.html")


@app.errorhandler(405)
def method_not_found(e):
    return render_template("405.html")


@app.route('/')
@login_required
def homepage():
    # user's directory
    userdir = "{}{}".format(Config.config['PATH']['output_dir'],
                            session['username'])
    model_list_file = "{}/modelo_lista.txt".format(userdir)
    vector_file = "{}/vectors.txt".format(userdir)
    paramfile = "{}/select_param.json".format(userdir)
    vectormapfile = "{}/vectorsmap.json".format(userdir)

    generalsolutionlist = '{}/resume.txt'.format(
        Config.config['PATH']['general_solution'])
    try:
        # if there's data load the model from the station requested by the user
        if os.path.isfile(model_list_file):
            lonlat, df = load_model(model_list_file)
            stations_param = np.array(df)

            # empty for the rest of the fields
            df_withvectors = copy.deepcopy(df)
            df_withvectors.assign(vector_e=np.nan, vector_n=np.nan,
                                  vector_z=np.nan, start_time_str=np.nan,
                                  end_time_str=np.nan)
            df_withvectors = np.array(df_withvectors)
            if os.path.isfile(vector_file):
                # load the vector data if there is a vector file
                middle = MiddleLayer(session['username'])
                df_withvectors = middle.middle_vectortable(df, vector_file)
        else:
            stations_param = []
            df_withvectors = []
            if os.path.isfile(generalsolutionlist):
                station_name = np.loadtxt(generalsolutionlist, usecols=[0],
                                          skiprows=1, dtype=str).tolist()
                lonlat = np.loadtxt(generalsolutionlist, usecols=[1, 2],
                                    skiprows=1).tolist()
                lonlat = [station_name, lonlat]
            else:
                lonlat = []
        # load the coordinates of the area in wich is going to be displayed
        # the vorticity field and the velocity vectors on the map
        if os.path.isfile(paramfile):
            with open(paramfile, 'r') as f:
                paramdata = json.load(f)
            try:
                fieldarea = [paramdata['Field longitude min'],
                             paramdata['Field longitude max'],
                             paramdata['Field latitude min'],
                             paramdata['Field latitude max']]
            except KeyError as e:
                fieldarea = []
            try:
                velocityarea = [paramdata['Vector longitude min'],
                                paramdata['Vector longitude max'],
                                paramdata['Vector latitude min'],
                                paramdata['Vector latitude max']]
            except KeyError as e:
                velocityarea = []
        else:
            fieldarea = []
            velocityarea = []

        if os.path.isfile(vectormapfile):
            with open(vectormapfile, 'r') as f:
                v2plot = json.load(f)
                v2plot = np.array([v2plot['x'], v2plot['y'], v2plot['vx'],
                                   v2plot['vy'], v2plot['ti'], v2plot['tf']]).T
                v2plot = v2plot.tolist()

        else:
            v2plot = []
        return render_template("index.html", stations=stations_param,
                               lonlat=lonlat, vectors=df_withvectors,
                               fieldarea=fieldarea, velocityarea=velocityarea,
                               vectors2plot=v2plot, databases=gnss_databases,
                               herenow=strftime("%Y-%m-%d %H:%M:%S", gmtime()))

    except TypeError as e:
        error = 'Fatal Exception'
        log = Logger()
        log.logger.error(PrintException())
        flash(error)
    return render_template("index.html", stations=[], lonlat=[], vectors=[],
                           fieldarea=[], velocityarea=[], vectors2plot=[],
                           herenow=strftime("%Y-%m-%d %H:%M:%S", gmtime()),
                           databases=gnss_databases)


@app.route('/select', methods=['POST'])
@login_required
def select():
    try:
        # manage the request of data from the user
        if request.method == 'POST' and request.form['btn_select'] == 'Select':
            latmin = float(request.form['latmin'])
            latmax = float(request.form['latmax'])
            lonmin = float(request.form['lonmin'])
            lonmax = float(request.form['lonmax'])
            ti = request.form['t_inicial']
            tf = request.form['t_final']

            middle = MiddleLayer(session['username'])
            middle.middle_select(lonmin, lonmax, latmin, latmax, ti, tf)
            flash('finished the data loading')
    except TypeError as e:
        error = 'Error. Please, check your query'
        log = Logger()
        log.logger.error(PrintException())
        flash(error)
    return redirect(url_for('homepage'))


@app.route('/about/')
def aboutpage():
    return render_template('about.html')


@app.route('/edit/<station>', methods=['POST', 'GET'])
@login_required
def edit(station):
    """
    Code for handle edit request
    """
    try:
        # transform input to a list with floats
        new_poly = int(request.form['poly'])
        new_jump = request.form['jump']
        new_fourier = request.form['fourier']
        new_logstart = request.form['logstart']
        new_logscale = request.form['logscale']

        middle = MiddleLayer(session['username'])
        middle.middle_edit(station, new_poly, new_jump, new_fourier,
                           new_logstart, new_logscale)

        flash("{} upgrated".format(station))
    except ValueError as e:
        error = 'Exception. Please, check your inputs (edit {}).\
If the error repeats contact the admin'.format(station)
        log = Logger()
        log.logger.error(PrintException())
        flash(error)

    return redirect(url_for('homepage'))


@app.route('/plot/<station>')
@login_required
def plots(station):
    """
    Plot request handlign
    """
    try:
        middle = MiddleLayer(session['username'])
        filename = middle.middle_plot(station)
        # cache_timeout should stop that t
        return send_file(filename)
    except (ValueError, IOError) as e:
        log = Logger()
        log.logger.error(PrintException())
        flash('Sorry, error ploting {}'.format(station))
        return 'error'


@app.route('/plot-v/<station>')
@login_required
def plots_v(station):
    """
    Plot request handlign
    """
    try:
        middle = MiddleLayer(session['username'])
        filename = middle.middle_plot(station, '-v')
        # cache_timeout should stop that t
        return send_file(filename)
    except (ValueError, IOError) as e:
        log = Logger()
        log.logger.error(e)
        flash('Sorry, error ploting {}'.format(station))
        return 'error'


@app.route('/vector/', methods=['POST', 'GET'])
@login_required
def calc_vector():
    """
    Calculate velocity vectors for the station 'station'
    """
    userdir = "{}{}".format(Config.config['PATH']['output_dir'],
                            session['username'])
    vector_file = "{}/vectors.txt".format(userdir)
    model_list_file = "{}/modelo_lista.txt".format(userdir)
    paramjsonfile = "{}/select_param.json".format(userdir)
    try:
        if request.method == 'POST' and request.form['btn_vector'] == 'Calculate':
            tv1 = request.form['t1']
            tv2 = request.form['t2']
            station = request.form['stationname']

            middle = MiddleLayer(session['username'])
            middle.middle_vector(station, tv1, tv2, vector_file,
                                 model_list_file, paramjsonfile)
            flash('{} vector sucessfull'.format(station))
    except (ValueError, TimeSeriesError) as e:
        if e is ValueError:
            log = Logger()
            log.logger.error(PrintException())
            flash('Start time should be lower than the end time')
        elif e is TimeSeriesError:
            log = Logger()
            log.logger.error(PrintException())
            flash('Time Series empty in the selected time')
    return redirect(url_for('homepage'))


@app.route('/download')
@login_required
def download():
    """
    zip the files in the user directory
    """
    output_dir = Config.config['PATH']['output_dir']
    os.chdir(output_dir)
    # check if zipfile exists
    if os.path.isfile('{}/yourcooldata.zip'.format(session['username'])):
        os.remove('{}/yourcooldata.zip'.format(session['username']))

    try:
        zf = zipfile.ZipFile('{}/yourcooldata.zip'.format(session['username']),
                             mode='w')
        zf.write('{}/series_lista.txt'.format(session['username']),
                 'list_series.txt')
        zf.write('{}/modelo_lista.txt'.format(session['username']),
                 'list_model.txt')

        for seriesfile in os.listdir('{}/series'.format(session['username'])):
            zf.write('{}/series/{}'.format(session['username'], seriesfile),
                     'series/{}'.format(seriesfile))
        for modelfile in os.listdir('{}/modelo'.format(session['username'])):
            zf.write('{}/modelo/{}'.format(session['username'], modelfile),
                     'model/{}'.format(modelfile))
        for pngfile in [png for png in os.listdir(session['username'])
                        if png.endswith('.png')]:
            zf.write('{}/{}'.format(session['username'], pngfile), pngfile)
        vectorfile = '{}/vectors.txt'.format(session['username'])
        if os.path.isfile(vectorfile):
            zf.write(vectorfile, vectorfile)
        zf.close()
        return send_file('{}{}/yourcooldata.zip'.format(output_dir,
                                                        session['username']),
                         as_attachment=True)

    except FileNotFoundError as e:
        log = Logger()
        log.logger.error(PrintException())
        flash('error retreaving the zip file')
        return redirect(url_for('homepage'))


@app.route('/field/',  methods=['POST'])
@login_required
def field():
    """
    Plot request handlign
    """
    userdir = "{}{}".format(Config.config['PATH']['output_dir'],
                            session['username'])
    vector_file = "{}/vectors.txt".format(userdir)
    model_list_file = "{}/modelo_lista.txt".format(userdir)
    if request.method == 'POST' and request.form['btn_field'] == 'Plot':
        latmin = float(request.form['latmin'])
        latmax = float(request.form['latmax'])
        lonmin = float(request.form['lonmin'])
        lonmax = float(request.form['lonmax'])
        ti = request.form['t_inicial']
        tf = request.form['t_final']

        grid_density = float(request.form['grid'])
        alfa = float(request.form['alfa'])
        # check input values
        if grid_density > 100000 or grid_density < 100:
            raise ValueError
        if latmin > latmax or lonmin > lonmax:
            flash('Lotitude or Latitude minimun value is bigger than the\
maximun value')
            raise ValueError
    middle = MiddleLayer(session['username'])
    middle.middle_field(model_list_file, vector_file, [lonmin, lonmax],
                        [latmin, latmax], [ti, tf], grid_density, alfa)

    return redirect(url_for('homepage'))



@app.route('/mapdata/<content>')
@login_required
def mapdata(content):
    output_dir = Config.config['PATH']['output_dir']
    if content == 'data':
        geojsonfile = "{}{}/out.geojson".format(output_dir,
                                                session['username'])
        if os.path.isfile(geojsonfile):
            mapdata = geojsonfile
            return send_file(mapdata, as_attachment=True)
        else:
            mapdata = ''
            return mapdata
    elif content == 'wz_figure':
        wzfigurefile = "{}{}/wz_field.png".format(output_dir,
                                                  session['username'])
        if os.path.isfile(wzfigurefile):
            mapdata = wzfigurefile
            return send_file(mapdata, as_attachment=True)
        else:
            mapdata = ''
            return mapdata

    elif content == 'wz_colorbar':
        colorbarfile = "{}{}/wz_colorbar.png".format(output_dir,
                                                     session['username'])
        if os.path.isfile(colorbarfile):
            mapdata = colorbarfile
            return send_file(mapdata, as_attachment=True)
        else:
            mapdata = ''
            return mapdata

    elif content == 'wk_figure':
        figurefile = "{}{}/wk_field.png".format(output_dir,
                                                session['username'])
        if os.path.isfile(figurefile):
            mapdata = figurefile
            return send_file(mapdata, as_attachment=True)
        else:
            mapdata = ''
            return mapdata

    elif content == 'wk_colorbar':
        colorbarfile = "{}{}/wk_colorbar.png".format(output_dir,
                                                     session['username'])
        if os.path.isfile(colorbarfile):
            mapdata = colorbarfile
            return send_file(mapdata, as_attachment=True)
        else:
            mapdata = ''
            return mapdata

    elif content == 's2_figure':
        figurefile = "{}{}/s2_field.png".format(output_dir,
                                                session['username'])
        if os.path.isfile(figurefile):
            mapdata = figurefile
            return send_file(mapdata, as_attachment=True)
        else:
            mapdata = ''
            return mapdata

    elif content == 's2_colorbar':
        colorbarfile = "{}{}/s2_colorbar.png".format(output_dir,
                                                     session['username'])
        if os.path.isfile(colorbarfile):
            mapdata = colorbarfile
            return send_file(mapdata, as_attachment=True)
        else:
            mapdata = ''
            return mapdata
    else:
        return ''


@app.route('/selectdb/', methods=['POST'])
@login_required
def select_db():
    """
    it makes the change between different GNSS databases
    """
    output_dir = Config.config['PATH']['output_dir']
    userdir = "{}{}".format(output_dir, session['username'])
    os.chdir(output_dir)

    try:
        if request.method == 'POST' and request.form['btn_db'] == 'Select DB':
            selected_db = request.form['select_db']
            Config.config['PATH']['general_solution'] = "{}general_solution_{}/".format(output_dir, selected_db)
            Config.config['PATH']['GPSdata'] = "txtfiles_{}".format(selected_db)
            Config.config['PATH']['ListaGPS'] = "station_list_{}.txt".format(selected_db)

            # save changes on ini file
            with open('./parameters.ini', 'w') as configfile:
                Config.config.write(configfile)
            if os.path.isdir(userdir):
                shutil.rmtree(userdir)

        return redirect(url_for('homepage'))

    except FileNotFoundError as e:
        log = Logger()
        log.logger.error(PrintException())
        flash('sorry, error selecting the DB')
        return redirect(url_for('homepage'))


# run function
def run():
    app.run()


if __name__ == "__main__":
    app.run()
