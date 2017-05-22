#! /usr/bin/python3

"""
@author: Vicente Yáñez
@date: 2017

Many auxiliary functions
"""

import json as json

import numpy as np


def dict2json(dictionary):
    json_obj = json.dumps(dictionary)
    return json_obj


def json2dict(json_str):
    dict_obj = json.loads(json_str)
    return dict_obj


def array2list(array):
    """
    Function that took one array and converts his binary string to
    a unicode string and save them in a new list.
    array input must be 1d array.

    Return 1d list
    """
    newlist = []
    for i in array:
        newlist.append(i.decode('utf-8'))

    return newlist
