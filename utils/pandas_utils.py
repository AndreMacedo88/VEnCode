#!/usr/bin/env python
# -*- coding: UTF-8 -*-

""" pandas_utils.py: utilities module for easier Pandas module integration """


def multi_set_data_frame(data, arrays, set_value):
    """
    Change multiple values on a pandas data frame in one line.
    :param data: pandas DataFrame. Note: must be of the same type as set_value.
            e.g.: DataFrame of "float" type and set_value=2.1.
    :param arrays: arrays of (row, col) directing to the cell to change the value.
    :param set_value: value to change to.
    """
    for i in arrays:
        data.at[i[0], i[1]] = set_value
