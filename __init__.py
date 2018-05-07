# flake8: noqa

import json
from os import path

fpath = path.abspath(path.dirname(__file__))

# initialize package version, content info, and revision uuid
with open(path.join(fpath, 'PKG_ID.json')) as vf:
    _pkginfo = json.load(vf)
__version__ = _pkginfo['version']
__content__ = _pkginfo['content']
__revisionid__ = _pkginfo['revisionid']

# verify external dependencies - raise an exception if required package is not installed
packages_to_test = ('pysb', 'numpy', 'netCDF4', 'messytables','pandas',
                     'xarray', 'dask', 'yaml', 'fastnumbers')

missing_dependencies = []
for dependency in packages_to_test:
    try:
        __import__(dependency)
    except ImportError as e:
        missing_dependencies.append(dependency)

if missing_dependencies:
    raise ImportError(
        "Missing required python modules(s){0}".format(missing_dependencies))

del dependency, missing_dependencies, packages_to_test, vf, fpath, json, path

# module level doc-string
__doc__ = """
BIS-OBIS - An OBIS processor library for Python,
A slice of the Biogeographic Information System
=====================================================================

This script downloads final source CSV files from the OBIS-USA collection, converts them to
netCDF and stages them for consumption by ERDDAP.  It also generates the datasets.xml ERDDAP
configuration file to serve those netCDF files as data sets.

Main Features
-------------
Functionality is limited (work in progress) to a few major tasks:

  - Retrieves source csv files and metadata information from ScienceBase
  - Processes source csv file using regex to determine formatting
  - Converts specified source CSV to NetCDF with a single "index" (row) dimension
  - Regenerates source CSV from NetCDF and compares csv files
  - Uses xarray-dask dataframe or messytables conversion approach
  - Controls memory consumption using data chunking
  - Incorporates basic error checking and message diagnostics
  - Evaluates dtypes, selects dtype using low-level logic + conflict resolution,
  - Summarizes results, rejects source file when encountering critical errors.
  - Outputs a basic processing report
  """
