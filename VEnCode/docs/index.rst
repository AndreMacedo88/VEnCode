Welcome to VEnCode's documentation!
===================================


.. toctree::
   :maxdepth: 3
   :caption: Contents:
   
   getting_started
   main_tools

Module for VEnCode-related projects based on FANTOM5 databases
--------------------------------------------------------------

|

.. image:: https://img.shields.io/pypi/v/VEnCode
    :target: https://pypi.org/project/VEnCode/
    :alt: PyPI
.. image:: https://img.shields.io/github/v/release/AndreMacedo88/VEnCode?include_prereleases
    :target: https://github.com/AndreMacedo88/VEnCode/releases
    :alt: GitHub release (latest by date including pre-releases)
.. image:: https://img.shields.io/github/release-date/AndreMacedo88/VEnCode
    :target: https://github.com/AndreMacedo88/VEnCode/releases
    :alt: GitHub Release Date
.. image:: https://travis-ci.com/AndreMacedo88/VEnCode.svg?branch=master
    :target: https://travis-ci.com/AndreMacedo88/VEnCode
.. image:: https://coveralls.io/repos/github/AndreMacedo88/VEnCode/badge.svg?branch=master
    :target: https://coveralls.io/github/AndreMacedo88/VEnCode?branch=master
.. image:: https://img.shields.io/pypi/pyversions/VEnCode
    :target: https://pypi.org/project/VEnCode/
    :alt: PyPI - Python Version
.. image:: https://img.shields.io/github/license/AndreMacedo88/VEnCode
    :target: https://github.com/AndreMacedo88/VEnCode/blob/Stable/LICENSE
    :alt: GitHub
.. image:: https://img.shields.io/github/issues/AndreMacedo88/VEnCode
    :target: https://github.com/AndreMacedo88/VEnCode/issues
    :alt: GitHub issues

|

This module contains classes and functions that perform intersectional genetics-related operations to find VEnCodes
using any matrix of cell types (columns) vs regulatory elements or markers (rows).

Moreover, it contains particular methods to make use of the databases provided by the `FANTOM5 consortium`_, namely the
CAGE enhancer and transcription start site (TSS) databases.

For more information on the VEnCode technology, please refer to `Macedo and Gontijo, GigaScience, 2020`_.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. Starting hyperlink targets:

.. _FANTOM5 consortium: http://fantom.gsc.riken.jp/5/data/
.. _Macedo and Gontijo, GigaScience, 2020: https://doi.org/10.1093/gigascience/giaa083