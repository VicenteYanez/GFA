#! /usr/bin/python3

import configparser
import os

"""
Class that loads the global parameters
"""


class Config():
    cwd = os.getcwd()
    paramfile = "{}/parameters.ini".format(cwd)

    config = configparser.ConfigParser()

    if os.path.isfile(paramfile):
        # user paramfile
        config.read(paramfile)
        # print("user file")
    else:
        # default config file
        thisfile = os.path.dirname(os.path.abspath(__file__))
        paramfile = "{}/parameters.ini".format(thisfile)
        config.read(paramfile)
        # print("default file")

    def main(self):
        print('nothing to do here...')


# test = Config()
# print(test.config['PATH']['eqfile'])
