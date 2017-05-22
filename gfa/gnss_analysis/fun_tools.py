#! /usr/bin/python
# -*- coding: utf-8 -*-

import re
import time
from datetime import datetime as dt


def toYearFraction(year, mes, dia, hora=0, minuto=0, segundo=0):
    """
    Funcion modificado de
    http://stackoverflow.com/questions/6451655/python-how-to-convert-datetime-dates-to-decimal-years
    Transforma fecha y hora a una fracción de un año.
    """

    def sinceEpoch(date):
        # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    date = dt(int(year), int(mes), int(dia))

    dia_fraction = (hora/24 + minuto/1440 + segundo/86400)/365.25

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction + dia_fraction


def string2date(fecha_texto):
    """
    Funcion que convierte una fecha en string a un objeto datetime
    """
    match = re.search(r'\d{2}-\d{2}-\d{4}', fecha_texto)
    fecha = dt.strptime(match.group(), '%d-%m-%Y').date()

    return fecha
