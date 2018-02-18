
"""
Class that manipulates the vector data loaded from vector.txt

Used a pandas dataframe
"""
from datetime import datetime, timedelta
from time import strftime
import copy
import json

import numpy as np
import pandas as pd


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


class VectorData():
    def __init__(self, vector_file):
        """
        Loads the vector file and join it with the rest of the station data
        """
        self.vecdf = pd.read_csv(vector_file)
        self.vector_file = vector_file

    def get_table(self, df):
        """
        Returns a numpy array with the data for the table showed in
        index.html
        """
        newdf = copy.deepcopy(self.vecdf)
        # change the format of date
        newdf['start_time'] = newdf['start_time'].apply(
            convert_partial_year)
        newdf['start_time_str'] = newdf['start_time'].apply(
            lambda x: datetime.strftime(x, '%Y-%m-%d'))
        newdf['end_time'] = newdf['end_time'].apply(
            convert_partial_year)
        newdf['end_time_str'] = newdf['end_time'].apply(
            lambda x: datetime.strftime(x, '%Y-%m-%d'))

        newdf['vector_n'] = newdf['vector_n'].apply(lambda x: '{0:.3f}'.format(
            float(x)))
        newdf['vector_e'] = newdf['vector_e'].apply(lambda x: '{0:.3f}'.format(
            float(x)))
        newdf['vector_z'] = newdf['vector_z'].apply(lambda x: '{0:.3f}'.format(
            float(x)))
        # join vectors vectors
        grouped = newdf.groupby('station')
        tempdf = grouped.aggregate({'vector_n': lambda x: tuple(x),
                                    'vector_e': lambda x: tuple(x),
                                    'vector_z': lambda x: tuple(x),
                                    'start_time_str': lambda x: tuple(x),
                                    'end_time_str': lambda x: tuple(x)})

        df_final = pd.merge(df, tempdf.reset_index(), on='station', how='left')
        df_final = df_final[['station', 'longitude', 'latitude', 'polynomial',
                             'jumps', 'fourier', 'log start', 'log scale',
                             'vector_e', 'vector_n', 'vector_z',
                             'start_time_str', 'end_time_str']]
        # format for output
        df_final['longitude'] = df_final['longitude'].apply(
            lambda x: '{0:.3f}'.format(float(x)))
        df_final['latitude'] = df_final['latitude'].apply(
            lambda x: '{0:.3f}'.format(float(x)))

        df_final['vector_e'] = df_final['vector_e'].apply(
            lambda x: str(x).replace(
                "(", "").replace(")", "").replace("'", "").replace(" ", ""))
        df_final['vector_n'] = df_final['vector_n'].apply(
            lambda x: str(x).replace(
                "(", "").replace(")", "").replace("'", "").replace(" ", ""))
        df_final['vector_z'] = df_final['vector_z'].apply(
            lambda x: str(x).replace(
                "(", "").replace(")", "").replace("'", "").replace(" ", ""))

        df_final['start_time_str'] = df_final['start_time_str'].apply(
            lambda x: str(x).replace(
                "(", "").replace(")", "").replace("'", "").replace(" ", ""))
        df_final['end_time_str'] = df_final['end_time_str'].apply(
            lambda x: str(x).replace(
                "(", "").replace(")", "").replace("'", "").replace(" ", ""))
        # format list result
        df_final = np.array(df_final)

        return df_final

    def select(self, df, lonrange, latrange, timerange):
        """
        returns arrays with vectors (positions, velocity and time) of the station
        selected by time and postion
        """
        newdf = copy.deepcopy(self.vecdf)
        newdf['start_time'] = newdf['start_time'].apply(
            convert_partial_year)
        newdf['end_time'] = newdf['end_time'].apply(
            convert_partial_year)
        df['longitude'] = df['longitude'].apply(
            lambda x: float(x))
        df['latitude'] = df['latitude'].apply(
            lambda x: float(x))
        selectdf = pd.merge(df, newdf, on='station', how='left')
        selectdf = selectdf[(selectdf['longitude'] > lonrange[0]) &
                            (selectdf['longitude'] < lonrange[1]) &
                            (selectdf['latitude'] > latrange[0]) &
                            (selectdf['latitude'] < latrange[1]) &
                            (selectdf['start_time'] > timerange[0]) &
                            (selectdf['end_time'] < timerange[1])]

        result1 = np.array(selectdf[['longitude', 'latitude', 'vector_e',
                                    'vector_n']]).T
        # result2 = np.array(selectdf[['start_time', 'end_time']]).T
        result21 = selectdf['start_time'].tolist()
        result22 = selectdf['end_time'].tolist()


        return result1[0], result1[1], result1[2], result1[3], result21, result22

    def check_time(self, station, ti, tf):
        """
        check if a time range is already inside the df
        return the new time range that is not in the dataframe
        return a empty list if there is a vector with the same time
        range for the station "station"
        if there is vectors in the old data that doesn't appear
        in the new list...they are returned in a separate list
        """
        list_start = self.vecdf[self.vecdf['station'] == station][
            'start_time'].tolist()
        list_end = self.vecdf[self.vecdf['station'] == station][
            'end_time'].tolist()

        old_times = np.array([list_start, list_end]).T
        print("old times: {}".format(old_times))
        input_times = np.array([ti, tf]).T
        print("input times {}".format(input_times))
        new_times = []
        remove_times = []
        # add the new input times to a list for calculation
        for times in input_times:
            if times not in old_times:
                new_times.append(times)
        # if a times in old_times not in the input data, add it
        # to a list for remove it later
        for times in old_times:
            if times not in input_times:
                remove_times.append(times)
        print("new times: {}".format(new_times))
        return new_times, remove_times

    def remove_vector(self, station, times2remove):
        """
        remove the vectors indicated
        """
        rmdf = copy.deepcopy(self.vecdf)
        for time2remove in times2remove:
            print("remove times: {}".format(time2remove))
            rmdf = rmdf[(rmdf['station'] != station) |
                        (rmdf['start_time'] != time2remove[0]) |
                        (rmdf['end_time'] != time2remove[1])]

        rmdf = rmdf[['station', 'vector_type', 'vector_e', 'vector_n',
                     'vector_z', 'c_e', 'c_n', 'c_z', 'start_time',
                     'end_time']]
        rmdf.to_csv(self.vector_file, index=False)
        return


def select_from_df(vector_file, lonrange, latrange, timerange):
    """
    returns arrays with vectors and positions of the station
    selected by time and postion
    """
    df = pd.read_csv(vector_file)
    df = df[(df['longitude'] > lonrange[0]) & (df['longitude'] < lonrange[1]) &
            (df['latitude'] > latrange[0]) & (df['latitude'] < latrange[1]) &
            (df['start_time'] > timerange[0]) & (df['end_time'] < timerange[1])
            ]
    result = np.array(df[['longitude', 'latitude',
                             'vector_e', 'vector_n']]).T

    return result[0], result[1], result[2], result[3]
