#! /usr/bin/python3

"""
@author: Vicente Yáñez

Simple script for config the log file
"""

import logging

# login config
logging.basicConfig(filename='../gfa.log', level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(name)s %(message)s')
logger = logging.getLogger(__name__)
