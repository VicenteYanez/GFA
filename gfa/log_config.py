#! /usr/bin/python3

"""
@author: Vicente Yáñez

Simple script for config the log file
"""

import logging

from gfa.load_param import Config


# login config
class Logger():
    def __init__(self):
        outputdir = Config.config['PATH']['output_dir']
        logging.basicConfig(filename='{}gfa.log'.format(outputdir),
                            level=logging.DEBUG,
                            format='%(asctime)s %(levelname)s %(name)s %(message)s')
        self.logger = logging.getLogger(__name__)
