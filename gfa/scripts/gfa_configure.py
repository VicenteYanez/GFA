
import configparser


"""
Function that creates a default parameters.ini file
in the current directory
"""


def main():
    config = configparser.ConfigParser()
    config['PATH'] = {'TimeSeries': './TimeSeries',
                      'ListaGPS': './TimeSeries/station_list.txt',
                      'GPSdata': './TimeSeries/txtfiles/',
                      'eqfile': './eq_file.txt',
                      'output_dir': './../output'}
    config['GNSS parameters'] = {'errmax': '10'}
    config['WEB'] = {'websession': 'False',
                     'offline_username': 'gfa'}
    with open('parameters.ini', 'w') as configfile:
        config.write(configfile)
