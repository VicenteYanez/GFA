
from datetime import datetime, timedelta
import time


def convert_partial_year(number):
    """
    Helper function to pass from fractional years to datetime object
    """
    year = int(number)
    d = timedelta(days=(number - year)*(365 + is_leap(year)))
    day_one = datetime(year, 1, 1)
    date = d + day_one
    return date


def is_leap(year):
    """
    Helper function than returns 1 if the year is leap, and 0 if it doesn't.

    it is use in convert_partial_year function
    """
    if year % 4 != 0:
        leap = 0
    elif year % 100 != 0:
        leap = 1
    elif year % 400 != 0:
        leap = 0
    else:
        leap = 1
    return leap


def fractlist2stryear(l):
    """
    Convert a list with year and fraction dates to a list with
    dates in string format %Y-%m-%d
    """
    stryear = [datetime.strftime(convert_partial_year(number), '%Y-%m-%d')
               for number in l]
    return stryear


def fractlist2dateobj(l):
    dateobj = [convert_partial_year(number) for number in l]
    return dateobj


def toYearFraction(date):
    def sinceEpoch(date):  # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction


def str2date(s):
    dateobj = datetime.strptime(s, "%Y-%m-%d")
    return dateobj


def str2yearfraction(s):
    yearfraction = toYearFraction(datetime.strptime(s, "%Y-%m-%d"))
    return yearfraction


def longstring2yearfraction(ls):
    """
    Transform a string with dates in the format %Y-%m-%d separated
    by a comma to a list with the dates in year and fraction
    """
    yrfr = [toYearFraction(datetime.strptime(s, "%Y-%m-%d"))
            for s in ls.split(',')]
    return yrfr


def list2yearfraction(l):
    """
    Transform a list with string dates in the format %Y-%m-%d to a list with
    the dates in year and fraction
    """
    yrfr = [toYearFraction(datetime.strptime(s, "%Y-%m-%d")) for s in l]
    return yrfr


def datelist2yearfraction(l):
    """
    Transform a list with string dates in the format %Y-%m-%d to a list with
    the dates in year and fraction
    """
    yrfr = [toYearFraction(dateobj) for dateobj in l]
    return yrfr


def datelist2strlist(l):
    """
    Transform a list with string dates in the format %Y-%m-%d to a list with
    the dates in year and fraction
    """
    strlist = [datetime.strftime(dateobj, '%Y-%m-%d') for dateobj in l]
    return strlist
