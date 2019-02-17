#!/usr/bin/env python
# -*- coding: UTF-8 -*-

""" general_util.py: General functions module """

import itertools as iter
import collections


def combinations_from_nested_lists(lst):
    """
    Generates tuples with combinations of one element of each list inside the first list.
    :param lst: the list to combine.
    """

    def helper(lst_):
        """
        Helps by allowing it's recursive use to generate a correctly shaped list to iter.product
        """
        if not any(isinstance(e, list) for e in lst_):
            lst_new.append(lst_)
        else:
            for z in lst_:
                if isinstance(z, list):
                    helper(z)
                else:
                    lst_new.append([z])

    if not any(isinstance(w, list) for w in lst):
        for g in lst:
            yield [g]
    else:
        lst_new = []
        for j in lst:
            if isinstance(j, list):
                helper(j)
            else:
                lst_new.append([j])
        lst = lst_new
        for i in iter.product(*lst):
            yield i


def flatten_irregular_nested_lists(l):
    """
    Flattens lists with sublists inside. sublists can be of multiple nesting levels.
    :param l: list to flatten
    """
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten_irregular_nested_lists(el)
        else:
            yield el


def e_value_normalizer(e_value_raw, k, n_celltypes):
    """
    Normalizes the e-value due to disparity in number of celltypes

    :param e_value_raw: value to normalize
    :param int k: number of rows, in practice it's the number of regulatory elements that give the VEnCode.
    :param int n_celltypes: Number of celltypes in the data (columns)
    :return: normalized e-value
    """
    coefs = {"a": -164054.1, "b": 0.9998811, "c": 0.000006088948, "d": 1.00051, "m": 0.9527, "e": -0.1131}
    e_value_expected = (coefs["m"] * k + coefs["e"]) * n_celltypes ** (
        coefs["d"] + ((coefs["a"] - coefs["d"]) / (1 + (k / coefs["c"]) ** coefs["b"])))
    e_value_norm = (e_value_raw / e_value_expected) * 100
    if e_value_norm < 100:
        return e_value_norm
    else:
        return 100
