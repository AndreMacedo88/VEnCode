#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
exceptions.py: handling custom exceptions.
"""


class NoVencodeError(Exception):
    """Raised when no VEnCode was found"""
    pass


class SampleTypeNotSupported(Exception):
    """Raised when sample type does ot exist in the FANTOM5 data"""
    def __init__(self, sample_type, cell_type, msg=None):
        if msg is None:
            msg = "sample_type - {} - currently not supported for cell type {}".format(sample_type, cell_type)
        super().__init__(msg)


class SampleTypeForCtpNotResolved(Exception):
    """Raised when sample type does ot exist in the FANTOM5 data"""
    def __init__(self, cell_type, msg=None):
        if msg is None:
            msg = "Not able to resolve sample_type for the cell type {}, please specify".format(cell_type)
        super().__init__(msg)
