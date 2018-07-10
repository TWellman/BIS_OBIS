
# coding: utf-8

# In[1]:

# %load obis_erddap_validate.py
# !/usr/local/bin/python3
# Framework: a processing component of the Biogeographic Information System (BIS)
# Group: Biogeographics Characterization Branch, USGS
# Authors: Collaboration (Tristan Wellman, BCB, USGS and John Long, USGS)
# Please see module information for additional details

"""

  This script downloads final source CSV files from the OBIS-USA collection, converts them to
  netCDF and stages them for consumption by ERDDAP.  It also generates the datasets.xml ERDDAP
  configuration file to serve those netCDF files as data sets.

  CSVs can be served from ERDDAP.  However, from the ERDDAP EDDTableFromAsciiFiles documentation:
      ASCII files are easy to work with, but they are not an efficient way to store/retreive data.
      For greater efficiency, save the files as NetCDF .nc files (with one dimension, "index",
      shared by all variables) instead.
"""
import logging
from pysb import SbSession
from datetime import datetime, timedelta, timezone
from xml.dom import minidom
import pkg_resources
import xml.etree.ElementTree  as ET
import json
import numpy
import netCDF4
import netCDF4 as netCDF
import messytables
import csv
import os
import sys
import glob
import json
import traceback
import decimal
import re
import getopt
import io
import pandas as pd
import xarray as xr
import collections
from itertools import zip_longest
from zipfile import ZipFile, is_zipfile  
from contextlib import ExitStack, redirect_stdout 
import dask
import yaml
from fastnumbers import fast_float as float_convert
from SPARQLWrapper import SPARQLWrapper, JSON
from urllib import request
try: 
    from BeautifulSoup import BeautifulSoup, Comment
except ImportError:
    from bs4 import BeautifulSoup, Comment


# In[2]:

# processing options ( ... most of them anyway)

def usage():
    print("""
    %s [options]
    

Options
-------------------

    --collectionid=<id> The ScienceBase Item ID of the collection to process.
    --sourcedir=<path> Directory in which to place source CSV fetch_csvs (defaults to ./source_data).
    --erddapdir=<path> Directory in which to place netCDF files (defaults to ./erddap_data/nc).
    --serverdir=<path> Directory in which netCDF files will be placed on the ERDDAP server.
      Defaults to /etddapData/nc.
    --recon_data_dir=<path> Directory in which to place reconstructed CSV created from NetCDF (defaults to ./erddap_data/nc/recreate_files).
    --report_dir=<path> Directory in which to place report comparison files (defaults to ./erddap_data/nc/processing_reports).
    --tempdir=<path> Directory in which to place temporary files (defaults to './erddap_data/nc/test/).
    --only_csv Download CSVs from ScienceBase (no netCDF or datasets.xml creation).
    --only_netcdf Create netCDF files from existing CSVs.
    --only_datasets_xml Create datasets.xml from existing netCDF files.
    --compare_csv2csv=<logical> whether to compare source csv to reconstructed csv generated from NetCDF (default = True) 
    --error_report=<logical> whether to output comparison report file when csv files are compared
    --postfilter=<logical> whether to post-filter content during file comparisons (default = False)
    --dump_csv=<logical> whether to delete regenerated csv file from netCDF, used for testing (default = True)
    --max_report=<integer> maximum errors to report per data column in processing reports (default = 20)
    --table_output=<logical> whether to print table summaries to screen (default = False) 
    
    --convert_method=<string> approach to process csv files("dataframe" or "messytables", default = dataframe)

    ** messytables method inputs **
    
    --virtual_datasets Create a virtual dataset for all netCDF files created from a single CSV.
    --window=<integer> Use the first n rows to guess CSV column datatypes (defaults to 500).
    --rowsperfile=<integer> Split the CSV data into netCDF files containing n rows each.
    --file_overwrite Overwrite existing source files.  Default is to skip files that exist locally.
    --proc_overwrite Overwrite existing converted files.  Default is to skip files that exist locally.
    --sample=<integer> Only process the first n rows of the CSV (for testing purposes).
    --verbose Additional logging.
    
    ** dataframe method inputs **
    
    --chunksize=<integer> reset maximum element chunk size (default = None)
    --prefilter=<logical> whether to prefilter source csv file, filters by "convert_chars" + regex functions
    --date_convert=<logical> whether to force format date variables (ex. use when ISO-8601 auto-parse fails, default = False) 
    --date_fmt=<string> string date format when force formatting is active (e.g. "%m-%d-%Y") 
    
    ** file comparison options **
    
    --int_float_accept=<logical> whether to allow integer-float equivalence, e.g. 1.00 = 1 (default = False) 
    
    
""" % (sys.argv[0]))


# In[3]:

#######################################################################
##                          program options                          ##
#######################################################################
    
# Retrieves default program options.
# The defaults may be overwritten using command line or *.yml file
# returns:  a dictionary of program options (commands)
# T. Wellman BCB Group, USGS and J. Long, USGS

def default_inputs():

    
    # ScienceBase Item ID - search for data files 
    #--------------------------------------------
    collection_id = '57fe93d5e4b0824b2d14cbe1' # '579b64c6e4b0589fa1c98118' #'57fe93d5e4b0824b2d14cbe1'  # '579b64c6e4b0589fa1c98118' 
    
    #
    # Dictionary of ScienceBase search terms - use list format only, not case sensitive 
    #
    file_srch = dict([('title', ['DarwinCore:event','DarwinCore:occurrence',
                                'DarwinCore:measurementOrFact','final processed source',
                                'Final Processed File']),
                      ('name', ['occurrence', 'event', 'measurementOrFact'])])
    file_srch['ftype_req'] = ['.csv', '.zip']


    #--------------------------------
    # Folder directories
    #--------------------------------

    #
    # main directory to run program (None defaults to current directory)
    #
    workdir = None

    #
    # folder with original (source) csv files
    #
    source_data_dir = './source_data'

    #
    # folder with netCDF files (converted from source csv)
    #
    erddap_data_dir = './erddap_data/nc_store'

    #
    # Relative path to the netCDF files on the destination server, for the datasets.xml file.
    #
    server_nc_directory = './erddap_data/nc_store'

    #
    # folder with regenerated csv files (converted from netCDF)
    #
    recon_data_dir = erddap_data_dir + '/recreate_files'

    #
    # folder to store error and other report files 
    #
    report_dir = erddap_data_dir + '/processing_reports'

    #
    # folder to store temporary files when data chunking active (dataframe option)
    #
    tempdir = erddap_data_dir + '/test/'
    
    #
    # path + name (vocabulary standard) in json file, "None" to bypass
    #
    vocab = pkg_resources.resource_filename(__name__, 'data/vocab_standards.json')
    vocab_name  = 'vocabulary standard'
    
    #--------------------------------
    # Processing method (options)
    #--------------------------------

    # Exploratory approaches at processing tabular data
    #  "dataframe" approach to process + convert csv files(xarray/dask/pandas) Tristan Wellman, USGS 
    #  "messytables" approach (csv_reader/messytables))  John Long, USGS
    #
    convert_method = "dataframe" # "messytables"  
    

    #--------------------------------
    # processing flags (options)
    #--------------------------------

    #
    # Whether to fetch metadata from ScienceBase.  If this is True and fetch_csvs
    # is false, it will only fetch metadata.
    #
    fetch_metadata = True 

    #
    # Whether to fetch source files from ScienceBase.
    #
    fetch_csvs = True

    #
    # Whether to create netCDF files from source csv files
    #
    create_netcdf_files = True

    #
    # Whether to create a datasets.xml from the netCDF files found in the erddap_data_dir
    #
    create_datasets_xml = True

    #
    # flag whether to overwrite source files.  If set to False, reprocessing input file will be skipped.
    #
    file_overwrite = False
    
    #
    # flag whether to overwrite converted files.  If set to False, reprocessing input file will be skipped.
    #
    proc_overwrite = False

    #
    # flag whether to regenerate test csv files from netCDF, compare original (source) csv to regenerated csv
    #
    compare_csv2csv = True

    #
    # flag whether to delete regenerated csv file from netCDF, after testing 
    #
    dump_csv = True

    #
    # flag whether to output comparison report file (if compare_csv2csv = True)
    #
    error_report = True

    #
    # flag to print comparison table summaries to screen
    #
    table_output = False

    #
    # maximum errors to report per data column in processing reports (one per file, if error_report active)
    #
    max_report = 50


    #--------------------------------
    # Messy table ONLY options
    #--------------------------------

    #
    # Whether to create a single virtual dataset from multiple netCDFs that were created from a single csv
    #
    create_virtual_datasets = False

    #
    # Whether to turn on verbose logging
    #
    verbose = False

    #
    # Size in number of rows of the sample.  Only takes effect if sample is greater than zero.
    #
    sample_size = 0

    #
    # Number of rows to use to guess column type
    #
    window = 500

    #
    # Max number of rows per netCDF file.  If greater than zero, multiple files will be created if necessary.
    # Zero Value signals to create a single netCDF file for the dataset.
    #
    rows_per_file = 0


    #--------------------------------
    # Dataframe ONLY options
    #--------------------------------
    
    #
    # force convert (date) variables to datetime objects 
    # 'datetime': iso datetime variable, 'string': iso string date, "integer": days since 1-1-1970 basedate, 
    # otherwise : do not convert)
    #
    date_convert = 'string' # 'datetime'
    
    #
    # date format when forced to change (date_convert = True)
    #
    date_fmt = '%Y-%m-%dT%H:%M:%SZ' # "%Y-%m-%d"  

    #
    # file processing chunk size (number of elements, "None" deactivates chunking) 
    #
    chunk_elements = 5e6
    
    #
    # attempt CF compliancy standards
    #
    cf_comply = False
    
    #
    # representation for missing string entries (if == '_FillValue', uses ' ' activates NetCDF _FillValue)
    #
    absent_string = 'NA'
    
    #
    #  specify netcdf type: NETCDF4, NETCDF4_CLASSIC, NETCDF3_64BIT, or NETCDF3_CLASSIC
    #
    netcdf_type = 'NETCDF4_CLASSIC'  
    
    #--------------------------------
    # file comparison options
    #--------------------------------

    #
    # file encoding format - incomplete application
    #
    string_fmt = 'UTF-8' # 'ISO-8859-1'
    
    #
    # flag whether to allow integer-float equivalence, e.g. 1.00 = 1 (yes) 1.0001 = 1 (no)
    #
    int_float_accept = False
    
    #
    # filename modification of source file, when source csv is reconcontructed from netCDF
    #
    fname_ext = '_redo_'
    
    
    #---------------------------------------
    # Misc. specifications (in progress)
    #---------------------------------------
    
    #
    # logging options (optional flag - log to screen, set log level (e.g. debug, info, warning)) 
    #
    log_screen = True
    log_level = logging.DEBUG
    
    #
    # used for testing only, limit processing to # datasets, default is false (off)
    proc_limit = False 

    # note: prefilter *currently* used in dataframe method, postfilter works for dataframe or messytable methods

    #
    # flag whether to prefilter source csv file, filters by "convert_chars" + regex functions 
    #
    prefilter = False

    #
    # flag whether to post-filter file reads during comparisons to allow differences, same technique as prefilter
    #
    postfilter = False
    
    
    #
    # qoute interpreter, default, if None --> autoconfigured
    #
    quote_style = csv.QUOTE_NONNUMERIC # csv.QUOTE_MINIMAL
    

    # dict of character strings to filter from dataset (prefilter or postfilter)
    # key, value : replaced term, modified term
    # note: search keys are not case sensitive 
    convert_chars =  collections.OrderedDict([
        ('\n' , ''),
        (',n/a,'  , ',' + absent_string + ','),
        (',none,'  , ',NA,'),
        (',na,'  , ',NA,') ])
    
    # Character string qoute format - used in Dataframe conversion method and file comparisons,
    # note: dataframe method re-evaluates during NetCDF processing, 
    # inputs below are defaults, may be updated or overwritten internally
    # 
    #
    q_fmt = ['"', '"{}"', "'"]
    qoute_format = '''[\,](?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)'''
    
    #
    # Dictionary of ERDDAP reserved variables (L.L.A.T) specifications
    #
    LLAT_specs = dict([
         ('longitude',{'destinationName':'longitude', 'units':'degrees_east'}),
         ('latitude',{'destinationName':'latitude', 'units':'degrees_north'}),
         ('altitude',{'destinationName':'altitude', 'units':'m', 'positive': 'up'}),
         ('depth', {'destinationName':'depth', 'units':'m', 'positive': 'down'}),
         ('eventDate',{'destinationName':'time'})  
     ])
    
    #
    # Search terms used on variable names, data type must be numeric or string 
    # (adhoc fuzzy logic, customize as needed)
    #
    numeric_terms = ['_meters', 'minimum_', 'maximum_', '_minimum_', '_maximum', 'decimal', 'in_meters']
    string_terms =  ['_id,', '_code']
    
    #
    # Dictionary of partial variable terms to infer units (adhoc logic, customize as needed)
    # key, value --> term, unit 
    variable_units = dict([
         ('inKg','kg'),
         ('in_meters','m'),
         ('inmeters','m'),
         ('WeightsN','N'),])
  
    # assemble commands dictionary - (arguments and options) 
    # ------------------------------------------------------

    commands = dict([('collection_id', collection_id),
    ('file_srch', file_srch), 
    ('workdir', workdir),
    ('source_data_dir', source_data_dir),
    ('erddap_data_dir', erddap_data_dir),
    ('server_nc_directory', server_nc_directory),
    ('recon_data_dir', recon_data_dir), 
    ('report_dir', report_dir), 
    ('tempdir', tempdir), 
    ('vocab', vocab),
    ('vocab_name', vocab_name),                
    ('convert_method', convert_method), 
    ('fetch_metadata', fetch_metadata), 
    ('fetch_csvs', fetch_csvs), 
    ('create_netcdf_files', create_netcdf_files),
    ('create_datasets_xml', create_datasets_xml), 
    ('file_overwrite', file_overwrite),
    ('proc_overwrite', proc_overwrite),                 
    ('compare_csv2csv', compare_csv2csv),
    ('dump_csv', dump_csv), 
    ('error_report', error_report), 
    ('table_output', table_output), 
    ('max_report', max_report), 
    ('create_virtual_datasets', create_virtual_datasets),
    ('verbose', verbose), 
    ('sample_size', sample_size), 
    ('window', window),
    ('rows_per_file', rows_per_file), 
    ('int_float_accept', int_float_accept), 
    ('date_convert', date_convert),
    ('numeric_terms', numeric_terms),                
    ('string_terms', string_terms),                 
    ('variable_units', variable_units),                 
    ('date_fmt', date_fmt),  
    ('string_fmt', string_fmt), 
    ('chunk_elements', chunk_elements),
    ('cf_comply', cf_comply),
    ('absent_string', absent_string), 
    ('netcdf_type', netcdf_type),                        
    ('fname_ext', fname_ext),
    ('log_screen', log_screen),
    ('log_level', log_level), 
    ('proc_limit', proc_limit),                     
    ('prefilter', prefilter), 
    ('postfilter', postfilter),
    ('quote_style', quote_style),
    ('convert_chars', convert_chars),
    ('q_fmt',q_fmt),   
    ('qoute_format', qoute_format),
    ('LLAT_specs' , LLAT_specs),
    ])
    
    return commands


# In[4]:


#######################################################################
#                       processing definitions 
####################################################################### 

#
# retrieve_source_files
#    sb:  ScienceBase pysb session
#    collection_id:  ScienceBase Item ID of the OBIS-USA source data collection
#    source_data_dir:  Directory in which to place downloaded source data
#    download:  Flag whether to download the file.  Download if True, otherwise only retrieve metadata
# Finds all CSVs marked as "Final Processed Source" in the given collection, retrieves metadata
# for them, and, optionally, downloads the file.
#
# Returns: Dict containing ScienceBase metadata keyed by source file name.
#
# John Long, developer USGS  
# T. Wellman, BCB Group, USGS
#

def retrieve_source_files(sb, commands):
    
    if commands['fetch_metadata'] == False and commands['fetch_csvs'] == False:
        
        # Evaluate local data files only, exit
        # ------------------------------------
        sources = {}
        listdir = [f for f in os.listdir(commands['source_data_dir']) if f.endswith('.csv')]
        for dfile in listdir:
            logger.debug("{}{}".format( dfile,' from different request'))
        return listdir, sources
            

    # read ScienceBase collection or single item
    # -------------------------------------------
    logger.info('** Searching ScienceBase records for data files and metdata **')
    
    results = sb.find_items({
        'ancestors': commands['collection_id'],
        'fields': 'title, body, files',
        'max': 10000})
    if not results['items']:
        single_item = sb.get_item(commands['collection_id'])
        results = {}
        results['items'] = [single_item]
        
    counter = dict([    
                    ('items to search', len(results['items'])),
                    ('total files', 0),
                    ('candidate files', 0),
                    ('downloaded files', 0),
                    ('scibase_size_mb', 0),
                    ('download_size_mb', 0),
                     ])

    # Retrieve files from ScienceBase 
    # -------------------------------
    source_files = {}
    no_item = []
    file_srch =  dict((k.lower(), [v.lower() for v in l ]  ) for k,l in commands['file_srch'].items())
    
    while results and 'items' in results:
        for item in results['items']:
            item_flag = False
            if 'files' in item:
                for item_file in item['files']:
                    counter['total files'] +=1
                    pass_flag = False                  
                    for key, terms in file_srch.items():
                        if key in item_file:
                            fname = item_file[key]
                            if fname:
                                fname = fname.lower()
                                if key == 'name':
                                    if any(t in fname for t in terms):
                                        if any(s in fname for s in file_srch['ftype_req']):
                                            pass_flag = True
                                elif any(t in fname for t in terms):
                                    pass_flag = True
                    if pass_flag:
                        item_flag = True
                        counter['scibase_size_mb'] += item_file['size']
                        counter['candidate files'] += 1
                        logger.debug('   {} {}'.format(counter['candidate files'], item_file['name']))
                        if '.zip' not in item_file['name']:
                            src_file = item_file['name']
                        else:
                            src_file = item_file['name'].split('.')[0] + '.csv'  
                        
                        if commands['fetch_csvs']:
                            sfile_path = os.path.join(commands['source_data_dir'], src_file)
                            if os.path.isfile(sfile_path):
                                tstamp = os.path.getmtime(sfile_path)
                                date_sfloc = datetime.fromtimestamp(tstamp)
                                try:
                                    date_sb = datetime.strptime(item_file['dateUploaded'],'%Y-%m-%dT%H:%M:%SZ')
                                except:
                                    try: 
                                        date_sb = datetime.strptime(item_file['dateUploaded'],'%Y-%m-%dT%H:%M:%S.%fZ')
                                    except:
                                        date_sb = False
                                if date_sfloc < date_sb - timedelta(hours=12):
                                    getflag = True
                                else:
                                    getflag = False
                                sbdate = date_sb.strftime('%Y-%m-%dT%H:%M:%SZ')
                                locdate = date_sfloc.strftime('%Y-%m-%dT%H:%M:%SZ')
                                logger.debug("    {}{}{}{}".format('SBase date: ', sbdate,', local date: ', locdate))
                            else:
                                getflag = True
                                
                            if not commands['file_overwrite'] and not getflag:
                                logger.debug("    {}".format('local source file is up to date'))
                                source_files[src_file] = item
                            else:
                                if '.csv' in item_file['name']:
                                    counter['downloaded files'] +=1
                                    counter['download_size_mb'] += item_file['size']
                                    sb.download_file(item_file['url'], item_file['name'],
                                                     commands['source_data_dir'])
                                    source_files[item_file['name']] = item
                                elif '.zip' in item_file['name']:
                                    zfile = os.path.join(commands['source_data_dir'], item_file['name'])
                                    if not os.path.isfile(zfile):
                                        sb.download_file(item_file['url'], item_file['name'],
                                                         commands['source_data_dir'])
                                    if is_zipfile(zfile):
                                        with ZipFile(zfile) as files_zipped:
                                            zip_objs = files_zipped.infolist()
                                            for zf in zip_objs:
                                                if any(s in zf.filename and 
                                                       s is not '.zip' for s in file_srch['ftype_req']):
                                                    counter['downloaded files'] +=1
                                                    dfile = os.path.join(commands['source_data_dir'], zf.filename)
                                                    try:
                                                        files_zipped.extract(member=zf, path=commands['source_data_dir']) 
                                                        counter['download_size_mb'] += os.path.getsize(dfile)
                                                        source_files[item_file['name']] = item
                                                        source_files[zf.filename] = source_files.pop(item_file['name'])
                                                    except:
                                                        no_item.append("\t{}{}".format('Unzip attempt failed for item: ', item['id']))
                                                        pass
                                                else:
                                                    logger.debug("{}{}".format(zf.filename, ' ignored - zipped file without .csv extension'))      
                                    else:
                                        logger.debug('warning - unzip attempt failed, structure was not recognized')
                                        no_item.append("\t{}{}".format('Unzip attempt failed for item: ', item['id']))
         
                        elif fetch_metadata:
                            if os.path.join(commands['source_data_dir'], src_file):
                                source_files[item_file['name']] = item
             
            if item_flag == False:
                no_item.append("\t{}{}".format('no applicable files found in item: ', item['id']))                                
        results = sb.next(results)
        
      
    # summarize results
    #------------------
    
    if no_item:
        logger.debug("{}".format('** Sciencebase items identified without candidate files **'))
        for i in no_item:
            logger.info(i)
        
    counter['download_size_mb'] = counter['download_size_mb']/1000000.0
    counter['scibase_size_mb']  = counter['scibase_size_mb']/1000000.0
    for k,v in counter.items():
        logger.info("   Sciencebase search summary: {}  {}".format(k,v))
        
    listdir = [f for f in os.listdir(commands['source_data_dir']) if f.endswith('.csv')]
    for source in source_files:
        if source not in listdir:
            logger.debug("{}{}{}".format('warning - ', source,' file from request is not in source file directory'))
    for dfile in listdir:
        if dfile not in source_files:
            logger.debug("{}{}".format(dfile,' from different request'))

    return listdir, source_files


# Purpose: filters out specified characters, whitespaces, double qoutes, etc
#          in data files row-by-row as single strings
# funtions:
#   a) uses regex operations for customizable character filtering
#   b) splits string using comma delimiters
# returns:
#   c) filtered list of comma delimited strings from single string
#
# T. Wellman, BCB Group, USGS

def filter_convert_row(textstr, func, regex_filter, regex_parse, match_indx):
    textstr =  regex_filter.sub(func, textstr)
    textstr =  regex_filter.sub(func, textstr)
    row_strp = [x.strip() for x in regex_parse.split(textstr)]
    row_list = [row_strp[r] for r in match_indx]
    return row_list


# Purpose: FOR CODE TESTING ONLY, brute force combine of multiple netcdf files
# funtions (caution - eager memory consumption!):
#   a) uses xarray operations for combine
#   b) loads as dask via chunking call
# returns:
#   c) combined file object
#
# T. Wellman, BCB Group, USGS
def combine_netcdfs(paths, chunk, dim, transform_func=None):
    def process_one_path(path):
        with xr.open_dataset(path) as ds:
            ds = ds.chunk(chunks=chunk, name_prefix='xarray-', token=None, lock=True)
            if transform_func is not None:
                ds = transform_func(ds)
            ds.load()
            return ds
    datasets = [process_one_path(p) for p in paths]
    combined = xr.concat(datasets, dim)
    return combined


# Purpose: Produces csv file from NetCDF using
# interpreted formatting, often derived from source data file,
# produced csv file can be used for basic QA/QC, specifically
# to evaluate data content in text form.
#
# # T. Wellman, BCB Group, USGS 
#
# Assumes csv file is well formed (within Pandas dataframe capabilities)
#
def regenerate_csv(fpaths, commands, err_msg = [], p_tasks = []):
    
    if not os.path.isfile(fpaths['rconpath']) or commands['proc_overwrite']: 

        # data processing step size - uses chunks or complete read
        # -------------------------------------------------------
        
        if 'chunksize' not in commands:
            if 'chunk_elements' in commands:
                commands['chunksize'] = max(int(commands['chunk_elements']/100), 1)
            else:
                commands['chunksize'] = False
        if commands['chunksize']:
            ds = xr.open_dataset(fpaths['ncpath'], chunks=commands['chunksize'], decode_times=False)
            ds = ds.chunk(chunks=commands['chunksize'], name_prefix='xarray-', token=None, lock=True)
        else:
            ds = xr.open_dataset(fpaths['ncpath'], chunks=commands['chunksize'], decode_times=False)
        for k in ds.dims:
            nrows = ds.dims[k]
            if nrows > 0:
                break
            else:
                return
        if commands['chunksize']:
            ndiv = int(numpy.ceil(nrows/commands['chunksize']))
            csize = commands['chunksize']
        else:
            ndiv = 1
            csize = nrows 
            

        # convert netcdf to csv using dataframe
        # ---------------------------------------
        try:
            os.remove(fpaths['rconpath'])
        except:
            pass
        
        if ndiv == 1:
            
            logger.debug("\t{}".format('Reconstructed source file from NetCDF in one read'))
            
            with open(fpaths['rconpath'], 'a') as csvfile:
                df = ds.to_dataframe()         
                df.columns = df.columns.str.strip()
                df.columns = df.columns.str.replace(r"[/]", '_per_')
                df.columns = df.columns.str.replace(r"[^A-Za-z0-9]+", '_')
                df.to_csv(csvfile, columns=list(df.columns.values), index=False,
                          quotechar=commands['q_fmt'][0],
                          quoting=commands['quote_style'],
                          encoding = commands['string_fmt'])
        else:
            
            logger.debug("\t{}{}".format('Reconstructed source file from NetCDF in chunks: ', ndiv))
            
            f = True
            with open(fpaths['rconpath'], 'a') as csvfile:
                for n in numpy.arange(0, ndiv):
                    if f:
                        dslice = ds.sel(index=slice(n*csize, (n+1)*csize-1))
                        df = dslice.to_dataframe()
                        df.columns = df.columns.str.strip()
                        df.columns = df.columns.str.replace(r"[/]", '_per_')
                        df.columns = df.columns.str.replace(r"[^A-Za-z0-9]+", '_')
                        df.to_csv(csvfile, columns=list(df.columns.values),
                                  encoding=commands['string_fmt'], index=False,
                                  quoting=commands['quote_style'], 
                                  quotechar=commands['q_fmt'][0])
                        f = False
                        
                    else:
                        dslice = ds.sel(index=slice(n*csize, (n+1)*csize-1))
                        df = dslice.to_dataframe()                                                    
                        df.to_csv(csvfile, header=False, index=False,
                                  encoding=commands['string_fmt'],
                                  quotechar=commands['q_fmt'][0],
                                  quoting=commands['quote_style'])
        ds.close()
        p_tasks.append('reconstructed csv')
        return

    

# Purpose: lazily compares two csv files line-by-line 
#          allows for different column ordering and number of rows,
#          assumes identical row ordering.
# functions:
#    a) attempts low-level logic to discern cause of disagreements
#    b) optionally performs character filtering and float-integer 
#        equivalance operations prior to data comparisons
#    c) summarizes file size and dimensional agreements 
#    d) summarizes entry-by-entry agreements 
#    e) summarizes error reductions after selected operations
#    f) optionally writes report file (summaries, individual errors)
#  produces:
#    g) table summaries of comparison diagnostics
#    h) optional report file

# T. Wellman, BCB Group, USGS

def compare_csv(fpaths, commands, err_msg = [], p_tasks = []):
    
    if not os.path.isfile(fpaths['cmp_report']) or commands['proc_overwrite']: 
        if os.path.isfile(fpaths['sfpath']) and os.path.isfile(fpaths['rconpath']):
            
            p_tasks.append('compared source csv to reconstructed csv')    
            fmt = "\n{}\n\n\t{}\n\t{}{}\t\n\t{}{}\n\n"
            hdr_msg = fmt.format('*** File comparison evaluation',
                                 'purpose: compare source and reconstucted *.csv files',
                                 'source file: ', fpaths['sfpath'], 'reconstructed file: ',
                                 fpaths['rconpath'])
            if commands['table_output']: logger.debug(hdr_msg) 
            f1_sz = int(os.path.getsize(fpaths['sfpath']))
            f2_sz = int(os.path.getsize(fpaths['rconpath']))
            if f1_sz == f2_sz:
                msg ='PASSED size comparison: file sizes are identical' 
            else:
                fmt = "\n\n{}\n\t{}{}\n\t{}{}\n\t{}{}"
                msg = fmt.format('** CAUTION ** file sizes are not identical ',
                                 'original file (bytes): ', f1_sz,
                                 'replicate file (bytes): ', f2_sz,
                                 'percent difference: ', float(numpy.divide(f2_sz - f1_sz, 0.01 * f1_sz)),2)
            err_msg.append(msg)

            # set-up regex functions to process text and use convert_chars dictionary
            # -----------------------------------------------------------------------
            regex_filter = re.compile("(%s)" % "|".join(map(re.escape, commands['convert_chars'].keys())), re.IGNORECASE) 
            func = lambda d: commands['convert_chars'][d.string[d.start():d.end()].lower()]    
            regex_parse = re.compile(commands['qoute_format'])


            # lazily process file content line-by-line - lower memory requirement
            # -------------------------------------------------------------------
            with ExitStack() as stack:
                
                file1 =  stack.enter_context(
                    open(fpaths['sfpath'],'r', encoding= commands['string_fmt'], errors='replace'))
                file2 =  stack.enter_context(
                    open(fpaths['rconpath'],'r', encoding= commands['string_fmt'], errors='replace')) 
                
                commands['QC_pass'] = True
                cntr = dict( [('o_all', 0),('r_all', 0),('share_all', 0),
                              ('no_fltr', 0),('no_fltr_typ', 0),('no_share', 0),
                              ('o_row', 0),('r_row', 0),('o_col', 0),('r_col', 0)])  

                # read column headers (assumes simple header labels, comma delimiters)
                # --------------------------------------------------------------------
                nf1 = next(file1)
                nf2 = next(file2)
                cntr['o_col'] = len(nf1.split(','))
                cntr['r_col'] = len(nf2.split(','))

                # extract header labels 
                splt = re.sub('''["']''', '', nf1).split(',')
                hdr_1 = [h.strip() for h in splt]
                columns_hdr = len(hdr_1)
                splt = re.sub('''["']''', '', nf2).split(',')
                hdr_2 = [h.strip() for h in splt]

                # find matching header labels, extra labels, and missing labels
                # -------------------------------------------------------------
                if hdr_1 == hdr_2: 
                    msg = '\n\n** AGREEMENT ** column labels match in sequence'
                    err_msg.append(msg)
                    absent = []
                    extra = []
                else:
                    msg = '\n\n** CAUTION ** column labels do not match in sequence, processing by matches'
                    err_msg.append(msg)
                    extra = list(set(hdr_2) - set(hdr_1))
                    absent = list(set(hdr_1) - set(hdr_2))
                    if extra:
                        msg = '\n\n{}'.format('** CRITICAL ERROR - extra data columns were identified' + str(extra))
                        err_msg.append(msg)
                        commands['QC_pass'] = False
                    if absent:
                        msg = '\n\n{}'.format('** CRITICAL ERROR - missing data columns were identified' + str(absent))
                        err_msg.append(msg)
                        commands['QC_pass'] = False
                    if not absent and not extra:
                        msg = '\n\t{}'.format('result - all data columns were matched, with different ordering')
                        err_msg.append(msg)

                match_f1    =  [hdr_1.index(f) for f in hdr_1 if f in hdr_2 ]
                match_label =  [hdr_1[f] for f in match_f1] 
                match_f2    =  [hdr_2.index(f) for f in match_label] 

                # read, process, and filter content, and diagnose issues line-by-line
                # -------------------------------------------------------------------
                errors = []
                fmt = '{}{}{}{}{}{}  {}'
                fmt_wc = "\n\t{}\n\t\t{}{}\n\t\t{}{}\n\t\t{}{}\n\t\t{}{}\n"
                err_counter = collections.defaultdict(int)
                for i in match_label:
                    err_counter[i] = 0

                if match_f1:
                    for row, (line1, line2) in enumerate(zip_longest(file1, file2)):
                        if line1 and line2:

                            # read raw data strings,
                            # ignore character string qoutes,
                            # --------------------------------------------
                            f1_splt = [x.strip() for x in regex_parse.split(line1)]
                            len_f1 = len(f1_splt)
                            f2_splt = [x.strip() for x in regex_parse.split(line2)]
                            len_f2 = len(f2_splt)

                            # check if overall number of entries match 
                            # -----------------------------------------
                            if len_f1 != columns_hdr or len_f2 != columns_hdr:
                                msg = '\t{}{}{}'.format('CRITICAL ERROR: ', 
                                                        ' row ', row,' columns do not match header columns')
                                err_msg.append(msg)
                                logger.critical(fmt_wc.format('CRITICAL ERROR identified when comparing files:',
                                      'data row number [count]: ', row,
                                      'data header columns [count]: ', columns_hdr,                                      
                                      'original file row columns [count]: ', len_f1,
                                      'reconstructed file row columns [count]: ', len_f2))
                                commands['QC_pass'] = False
                                return
                            
                            # reduce comparisons to matching columns
                            # ---------------------------------------
                            f1_raw = [f1_splt[f] for f in match_f1]
                            f2_raw = [f2_splt[f] for f in match_f2]
                            
                            # optional - post-filter data prior to comparisons
                            # ------------------------------------------------
                            if commands['postfilter']:
                                f1 = filter_convert_row(line1, func, regex_filter, regex_parse, 
                                                        match_f1)
                                f2 = filter_convert_row(line2, func, regex_filter, regex_parse, 
                                                        match_f2)

                                # evaluate agreement by row (string) of data
                                # ------------------------------------------
                                if commands['int_float_accept']:
                                    for j,f in enumerate(f1_raw):
                                        cntr['share_all'] += 1
                                        if f1_raw[j] != f2_raw[j]:
                                            cntr['no_share'] += 1
                                        if f1[j] != f2[j]:
                                            cntr['no_fltr'] += 1
                                        if float_convert(f1[j]) != float_convert(f2[j]):    
                                            cntr['no_fltr_typ'] += 1
                                            err_counter[match_label[j]] +=1
                                            if err_counter[match_label[j]] <= commands['max_report']:
                                                err = fmt.format('Error row: ', row, ' col: ', match_label[j],
                                                             '  original : reproduced: ', str(f1_raw[j]), str(f2_raw[j]))
                                                errors.append(err)
                                else:
                                    for j,f in enumerate(f1_raw):
                                        cntr['share_all'] += 1
                                        if f1_raw[j] != f2_raw[j]:
                                            cntr['no_share'] += 1
                                        if f1[j] != f2[j]:
                                            cntr['no_fltr'] += 1
                                            cntr['no_fltr_typ'] += 1
                                            err_counter[match_label[j]] +=1
                                            if err_counter[match_label[j]] <= commands['max_report']:
                                                err = fmt.format('Error row: ', row, ' col: ', match_label[j],
                                                             '  original : reproduced: ', str(f1_raw[j]), str(f2_raw[j]))
                                                errors.append(err)
                            else:
                                if commands['int_float_accept']:
                                    for j,f in enumerate(f1_raw):
                                        cntr['share_all'] += 1
                                        if f1_raw[j] != f2_raw[j]:
                                            cntr['no_share'] += 1
                                            cntr['no_fltr'] += 1
                                            if float_convert(f1_raw[j]) != float_convert(f2_raw[j]):    
                                                cntr['no_fltr_typ'] += 1
                                                err_counter[match_label[j]] += 1
                                                if err_counter[match_label[j]] <= commands['max_report']:
                                                    err = fmt.format('Error row: ', row, ' col: ', match_label[j],
                                                                 '  original : reproduced: ', str(f1_raw[j]), str(f2_raw[j]))
                                                    errors.append(err)
                                else:
                                    for j,f in enumerate(f1_raw):
                                        cntr['share_all'] += 1
                                        if f1_raw[j] != f2_raw[j]:
                                            cntr['no_share'] += 1
                                            cntr['no_fltr_typ'] += 1
                                            cntr['no_fltr'] += 1
                                            err_counter[match_label[j]] += 1
                                            if err_counter[match_label[j]] <= commands['max_report']:
                                                err = fmt.format('Error row: ', row, ' col: ', match_label[j],
                                                             '  original : reproduced: ', str(f1_raw[j]), str(f2_raw[j]))
                                                errors.append(err)
                            cntr['o_row'] += 1
                            cntr['r_row'] += 1
                        else:
                            if line1 and not line2:
                                cntr['o_row'] += 1
                            else:
                                cntr['r_row'] += 1  

                    cntr['o_all'] = numpy.multiply(cntr['o_row'], cntr['o_col'])
                    cntr['r_all'] = numpy.multiply(cntr['r_row'], cntr['r_col'])


                    # output summary tables to screen and/or report file
                    # ---------------------------------------------------
                    if commands['QC_pass'] and (commands['table_output'] or commands['error_report']):

                        # table comparing overall file dimensions and column labels
                        # ---------------------------------------------------------
                        upper_label  = "\n\n{}".format('Summary table: overall dataset dimensions comparison') 
                        output_format = ['rows', 'columns']
                        measures = ['source file', 'reconstructed file', '(row, col) difference']
                        data = numpy.array([[cntr['o_row'], cntr['o_col']],
                                        [cntr['r_row'], cntr['r_col']],
                                        [cntr['o_row']-cntr['r_row'],
                                         cntr['o_col']-cntr['r_col']]])
                        overall_cmp = pd.DataFrame(data, measures, output_format)

                        # simple table summarizing comparison of individual data entries
                        # -------------------------------------------------------
                        middle_label = "\n\n{}".format('Summary table: entry-by-entry data agreement')
                        entry_measures = ['total source entries', 'shared (common) entries', 'matching of shared entries', 
                                          'non-matching of shared entries', 'missing entries', 'extra entries']
                        output_format = ['count', 'percentage']
                        absent_data = numpy.multiply(len(absent), cntr['o_row'])
                        extra_data  = numpy.multiply(len(extra), cntr['r_row']) 
                        share_match = cntr['share_all'] - cntr['no_share']
                        data = numpy.array([[cntr['o_all'], 'N/A'],
                                   [cntr['share_all'] , numpy.divide( cntr['share_all'], cntr['o_all']) * 100.0],
                                   [share_match, numpy.divide(share_match,cntr['share_all']) * 100.0],
                                   [cntr['no_share'], numpy.divide(cntr['no_share'],cntr['share_all']) * 100.0],
                                   [absent_data, numpy.divide(absent_data, cntr['o_all']) * 100.0 ],         
                                   [extra_data, numpy.divide(extra_data, cntr['r_all']) * 100.0 ],
                                   ])
                        entry_cmp = pd.DataFrame(data, entry_measures, output_format)

                        # simple table of comparisons after character filtering and equivalence operations  
                        # -------------------------------------------------------------------------
                        if cntr['no_share'] > 0:
                            msg = '\n\n{}'.format('** CAUTION ** found entry-by-entry disagreements (non-matches):')
                            err_msg.append(msg)
                            if cntr['no_fltr_typ'] == 0:
                                msg = '\n\n{}'.format('** PASSED ** after filtering and equivalence removed errors')
                                err_msg.append(msg)
                            else: 
                                msg = '\n\n{}'.format('** CAUTION ** errors remain after post-filter and equivalence allowances')
                                err_msg.append(msg)  
                            data = numpy.array([
                                [cntr['no_share'], numpy.divide(cntr['no_share'], cntr['share_all']) * 100.0],
                                [cntr['no_fltr'], numpy.divide(cntr['no_fltr'], cntr['share_all']) * 100.0],
                                [cntr['no_fltr_typ'], numpy.divide(cntr['no_fltr_typ'], cntr['share_all']) * 100.0],]) 
                        else:

                            # simple tale of error reductions after filtering and allowing float integer equivalence
                            # ----------------------------------------------------------------------------
                            data = numpy.array([[cntr['no_share'], 'N/A'],
                                [cntr['no_fltr'], 'N/A'],[cntr['no_fltr_typ'], 'N/A'],])   
                            msg = '\n\n{}'.format('** PASSED ** entry-by-entry values in agreement after operations')
                            err_msg.append(msg)

                        lower_label  = "\n\n{}".format('Summary table: discrepancy reductions from post-filtering and' + 
                                                       ' float-integer equivalence allowances')
                        error_measures = ['non-matching', 'filtering', 'filtering + equivalence']   
                        output_format = ['count', 'discrepancies remaining [%]']     
                        error_cmp = pd.DataFrame(data, error_measures, output_format)


                    # simple print messages/reports table to screen - optional 
                    # -------------------------------------------
                    if commands['table_output'] and commands['QC_pass']:
                        for e in err_msg:
                            logger.debug(e)  
                        logger.debug(upper_label); display(overall_cmp)
                        logger.debug(middle_label); display(entry_cmp)
                        logger.debug(lower_label); display(error_cmp)


                    # simple error report file of entry by entry conflicts - optional
                    # --------------------------------------------------------
                    if commands['error_report']:
                        
                        with open(fpaths['cmp_report'], "w") as ewrite:
                            disclaimer = "{}\n\t{}\n\t{}\n\t{}\n\t{}\n\t".format(
                                '*** Purpose: This file provides a basic summary of results while processing a *.csv data file, identified below.',
                                'The *.csv data file was evaluated for content and converted into NetCDF, considered viable for data release.', 
                                'The NetCDF file was then converted back to *.csv format, in most cases with minor adjustments to encoding and',
                                'data content. Validation comparisons were performed between the original source csv and regenerated csv files.',
                                'Results of the comparisons and other relevant information are shown below in preliminary form.')
                            agency = "{}\n\t{}".format('*** Processing agency: Biogeographic Characterization Branch,', 
                                                             'Core Science Systems, U.S. Geological Survey')
                            underscore = '--------------------------------------------------------------'
                            date_stamp = '*** Date-time of operation: ' + str(datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3])
                                
                            if commands['QC_pass']:
                                msg = '\n\n** TASKS COMPLETED:'
                                p_tasks.append('created report file')
                                for p in p_tasks:
                                    msg = msg + str(("\n\t{}".format(p)))
                                err_msg.append(msg + '\n\n')
                                entry_header = 'Unresolved comparisons AFTER optional post-filter/equivalence operations'
                                fmt = "{}\n\n{}\n\n{}\n{}\n\n{}\n{}\n{}\n\n\n{}\n{}\n{}\n\n\n{}\n{}\n{}\n\n\n"
                                ewrite.write(fmt.format(disclaimer, agency, hdr_msg, date_stamp, upper_label, underscore, 
                                                    overall_cmp.to_string(), middle_label, underscore, entry_cmp.to_string(),
                                                    lower_label, underscore, error_cmp.to_string(),))
                                ewrite.write("{}\n{}\n".format('Messages:', underscore))
                                for e in err_msg:
                                    ewrite.write(e)
                                ewrite.write("\n\n{}\n{}\n".format(entry_header, underscore))
                                ewrite.write("\n{}{}{}\n\n".format('output of individual errors limited to <=',
                                                                      commands['max_report'], ' per column'))
                                for err in errors:
                                    ewrite.write(err + '\n')
                                    
                            else:
                                fmt = "{}\n\n{}\n{}\n\n"
                                ewrite.write(fmt.format(disclaimer, hdr_msg, date_stamp))
                                underscore = '--------------------------------------------------------------'
                                ewrite.write("{}\n{}\n".format('Messages:', underscore))
                                for e in err_msg:
                                    ewrite.write(e)              
                    else:
                        logger.debug("\t{}{}".format('Print error report is turned off, error_report = ', commands['error_report']))
                else:           
                    logger.critical('Error - no matching column labels between datsets - terminating routine')
                    return
                
            if commands['dump_csv']:
                if os.path.exists(fpaths['rconpath']): 
                    logger.debug("\t{}".format('removed reconstructed source file'))
                    os.remove(fpaths['rconpath']) 
            return
        else:
            logger.warning("\t{}".format('warning **aborting comparison** files not found'))
            commands['QC_pass'] = False
            
        if commands['dump_csv']:
            if os.path.exists(fpaths['rconpath']): 
                logger.debug("\t{}".format('removed reconstructed source file'))
                os.remove(fpaths['rconpath']) 

    commands['QC_pass'] = True
    return

#
# csv_to_nc manual method
#     file_name: CSV (or ZIP containing a CSV) file to parse
#     source_data_dir: Directory containing CSV source files
#     metadata: ScienceBase item metadata for the file
#     erddap_data_dir: Directory in which to write the output netCDF file
#     rows_per_file
#
#     John (Dell) Long conversion approach, USGS
#
# Convert the given CSV to a netCDF with a "row" dimension.
#
def csv_to_nc_messy(file_name, source_data_dir, metadata, commands):
    
    logger.info("\t{}".format('Implementing messytables file conversion method'))
    commands['QC_pass'] = True # currently QC pass forced as True, functionality not available for messytables

    base_name, extension = os.path.splitext(file_name)
    full_file_name = os.path.join(source_data_dir, file_name)
    if 'csv' in extension.lower() or 'zip' in extension.lower():
        fpath = os.path.join(commands['erddap_data_dir'], base_name + '*.nc')
        file_exists = len(glob.glob(os.path.join(commands['erddap_data_dir'], base_name + '*.nc'))) > 0

        if not commands['proc_overwrite'] and file_exists:
            print('Not writing netCDF for %s, already exists' % (full_file_name))
            return
        (row_set, types, fh) = process_csv(full_file_name, extension, commands['window'])
        try:
            (variables, variable_data, num_rows) = normalize_data_for_netcdf(row_set, types, commands)
            offset = 0
            if commands['rows_per_file'] > 0:
                page_size = commands['rows_per_file']
            else:
                page_size = num_rows
            while offset < num_rows:
                create_netcdf(commands, base_name, metadata, variables, variable_data, offset, page_size, num_rows)
                offset += page_size
        except:
            raise
        finally:
            if fh:
                fh.close()
    else:
        logger.info('Unknown file %s type (%s), skipping' % (full_file_name, extension))

#
# process_csv
#     full_file_name: Path and filename of the csv
#     extension: File extension.  Can be "zip"" or "csv."
#
# Process the CSV, guessing column data types.
#
# Returns: The proceseed row_set, a dictionary of column types keyed by column name, and the
# open file handle (use this to close the file once the row_set is processed).
#
def process_csv(full_file_name, extension, window):
    row_set = []
    types = []
    fh = None
    logger.debug('\tProcessing %s' % (full_file_name))
    try:
        fh = open(full_file_name, 'rb')
        # Parse the CSV into a table set.  CSVs only have one table, so grab the first.
        if 'csv' in extension.lower():
            table_set = messytables.CSVTableSet(fh, window=window)
        elif 'zip' in extension.lower():
            table_set =  messytables.zip.ZIPTableSet(fh, window=window)

        row_set = table_set.tables[0]

        # Grab the header information, and set the iterator past it.
        offset, headers = messytables.headers_guess(row_set.sample)
        row_set.register_processor(messytables.headers_processor(headers))
        row_set.register_processor(messytables.offset_processor(offset + 1))

        # Guess column types and tell the row set to apply these types
        types = messytables.type_guess(row_set.sample, strict=True)
        row_set.register_processor(messytables.types_processor(types, strict=False))        
    except:        
        traceback.print_exc(file=sys.stdout)
        logger.warning('An error occurred, skipping %s' % (full_file_name))
        if fh:
            fh.close()
        row_set = types = []
        fh = None
    return (row_set, types, fh)
    

#
# normalize_data_for_netcdf
#     row_set: Messytables Rowset to process_csv
#     types: Dictionary of column types keyed by column NameError
#
# Processes the Messytables Rowset, normalizing the values in each cell for the column data type.
#
# Returns: List of variables containing interesting data, a dictionary of variable data keyed bytes
# column type, and the number of rows in the CSV.
#
def normalize_data_for_netcdf(row_set, types, commands):

    # Iterate the data and store for pre-processing
    variables = []
    variable_data = {}
    num_rows = 0

    for row in row_set:
        num_rows += 1
        for column_index, cell in enumerate(row):
            if cell.column not in variable_data:
                # Special case:  ERDDAP requires that altitude and depth both be numeric
                if cell.column.lower() in ['altitude', 'depth']:
                    types[column_index] = messytables.types.DecimalType()
                variables.append((cell.column, types[column_index]))
                variable_data[cell.column] = []
            variable_data[cell.column].append(nc_value(cell, types[column_index], commands))
        if commands['sample_size'] > 0 and num_rows == commands['sample_size']:
            break

    interesting_variables = []
    for variable, variable_type in variables:
        if has_values(variable_type, variable_data[variable]):
            if isinstance(variable_type, messytables.types.StringType):
                strlen = len(max(variable_data[variable], key=len))
                new_data = []
                for val in variable_data[variable]:
                    # Strip out non-ascii characters
                    ve = val.encode('ascii', errors='ignore').decode('ascii')
                    # Convert string to fixed length character array
                    new_data.append(netCDF4.stringtoarr(ve, strlen))
                variable_data[variable] = new_data
            interesting_variables.append((variable, variable_type))

    return (interesting_variables, variable_data, num_rows)

#
# create_netcdf
#     erddap_data_dir: Directory in which to write the netCDF files
#     base_name: File base name, used as the dataset name
#     metadata: ScienceBase metadata for the dataset
#     variables: List of tuples containing column names and types
#     variable_data: Dictionary of processed variable data, keyed by column name
#     offset: Index of starting row
#     page_size: Number of rows to write to this netCDF file
#     num_rows: Total number of rows in the data
#
# Creates a netCDF file from the processed CSV data.  Accepts an offset and page_size to allow multiple
# files to be written for each dataset.
#
def create_netcdf(commands, base_name, metadata, variables, variable_data, offset, page_size, num_rows):

    if num_rows > page_size:
        base_name = "%s_part%d" % (base_name, offset / page_size + 1)
    nc_fname = os.path.join(commands['erddap_data_dir'], base_name + ".nc")
    end_row = min(offset + page_size, num_rows)
    logger.debug("\tCreating %s for rows %d to %d" % (nc_fname, offset + 1, end_row))

    try:
        rootgrp = netCDF4.Dataset(nc_fname, "w", format="NETCDF3_CLASSIC")
        rootgrp.set_fill_on()

        # Set metadata attributes
        if metadata:
            rootgrp.setncatts(metadata)
            
        index_name = 'index'# row dimension name
            
        # Create netCDF dimensions for number of rows, and the size of each string field
        # row = rootgrp.createDimension('row', end_row - offset)
        row = rootgrp.createDimension(index_name, end_row - offset)  # Tristan change
        i = 0
        dims = {}
        uninteresting_vars = []

        for variable, variable_type in variables:
            if isinstance(variable_type, messytables.types.StringType):
                # Create a dimension for the character array length
                i += 1
                dim_name = "STRING%d" % (i)
                dims[variable] = dim_name
                # Use max length of all data so sizes stay constant across multiple netCDF files
                rootgrp.createDimension(dim_name, max([len(max(variable_data[variable], key=len)), 1]))

        # Create netCDF variables for each column.  We have to do this in a second pass, because all dimensions must
        # exist before we create variables.
        for variable, variable_type in variables:
            variable_data_slice = variable_data[variable][offset:end_row]
            if commands['verbose']:
                logger.debug('variable(%s)(%s)' % (variable, variable_type))
            if isinstance(variable_type, messytables.types.StringType):
                nc_var = rootgrp.createVariable(variable,nc_type(variable_type),(index_name,dims[variable],))
                nc_var[:] = numpy.asarray(variable_data_slice)
            else:
                nc_var = rootgrp.createVariable(variable,nc_type(variable_type),(index_name,))
                nc_var[:] = variable_data_slice
        rootgrp.close()
    except:
        logger.warning("Error creating netCDF %s, removing" % (nc_fname))
        traceback.print_exc(file=sys.stdout)
        os.remove(nc_fname)
#
# has_values
#     variable_type: Variable type
#     variable_data: Variable data
#
# Return whether the variable data contains data besides fill vaules
# Does not address data content
def has_values(variable_type, variable_data):
    ret = False
    if not isinstance(variable_type, messytables.types.BoolType):
        for var_value in variable_data:
            if var_value is not netCDF4.default_fillvals[nc_type(variable_type)]:
                ret = True
                break
    return ret

#
# nc_value
#     cell:  Messytables cell
#     column_type:  Type of the column the cell is in
#
# Return the given value in a format suitable for netCDF4
# TODO: Add nc attribute to netCDF for date mask.
#
def nc_value(cell, column_type, commands):
    ret = netCDF4.default_fillvals[nc_type(column_type)]
    try:
        if isinstance(column_type, messytables.types.DateType):
            ret = cell.value.date().isoformat()
        elif isinstance(column_type, messytables.types.StringType):
            if cell.value not in ['NA']:
                ret = str(cell.value).strip()
            else:
                ret = ''
        elif isinstance(column_type, messytables.types.IntegerType):
            ret = int(float(cell.value))
        elif isinstance(column_type, messytables.types.DecimalType):
            ret = float(cell.value)
        elif isinstance(column_type, messytables.types.BoolType):
            ret = str(cell.value)
        elif cell.value:
            ret = str(cell.value)
    except:
        if commands['verbose'] and cell.value is not None:
            logger.debug('filling %s(%s) column with %s for %s' % (cell.column, str(column_type), str(ret), str(cell.value)))
    finally:
        return ret

#
# nc_type
#     orig_type:  Type from messytables to convert to netCDF type
#
# Returns the netCDF type corresponding to the given messytables type
#
def nc_type(orig_type):
    ret = 'S1'
    if isinstance(orig_type, messytables.types.IntegerType):
        ret = 'i4'
    elif isinstance(orig_type, messytables.types.DecimalType):
        ret = 'f8'
    elif isinstance(orig_type, messytables.types.DateType):
        ret = 'S1'
    elif isinstance(orig_type, messytables.types.StringType):
        ret = 'S1'
    elif isinstance(orig_type, messytables.types.BoolType):
        ret = 'S1'
    return ret
    
#
# print_object
#     o: Object to print
# Prints the input object as comma-separated values (useful for debugging)
#
def print_object(o):
    logger.debug(', '.join("%s: %s" % item for item in vars(o).items()))

#
# write_datasets_xml
#     erddap_data_dir: Data containing the netCDF files
#     config_dir: Directory to which to write the datasets.xml file
#     create_virtual_datasets: Whether to combine multiple netCDF files for the
#         same dataset into one virtual ERDDAP dataset.
#
# Write the datasets.xml file for the given netCDF files.
#
def write_datasets_xml(commands, config_dir):
    logger.info('{}'.format('** Writing Datasets.xml file for ERDDAP'))
    datasets_xml_root = create_datasets_xml_root()
    processed = []
    r = re.compile(r'^(.*)_part\d+$')
    for file_name in glob.glob(os.path.join(commands['erddap_data_dir'], '*.nc')):
        base_name, extension = os.path.splitext(os.path.basename(file_name))
        m = re.match(r, base_name)
        if commands['create_virtual_datasets'] and m:
            dataset_name = m.group(1)
        else:
            dataset_name = base_name
        if dataset_name not in processed:
            logger.debug('{}{}'.format('   writing in dataset ', dataset_name))
            add_dataset(datasets_xml_root, dataset_name, file_name, commands['server_nc_directory'], 
                        commands['create_virtual_datasets'])
            processed.append(dataset_name)

    xmlstr = minidom.parseString(ET.tostring(datasets_xml_root)).toprettyxml(indent="   ", encoding="UTF-8")
    with open(os.path.join(config_dir, "datasets.xml") , "wb") as f:
        f.write(xmlstr)
                
#
# create_datasets_xml_root
#
# Create the DOM and add top level elements for the datasets.xml file
#
def create_datasets_xml_root():
    root = ET.Element('erddapDatasets')

    comment = ET.Comment('Generated from OBIS sources in ScienceBase')
    root.append(comment)

    convertToPublicSourceUrl = ET.SubElement(root, 'convertToPublicSourceUrl')

    requestBlacklist = ET.SubElement(root, 'requestBlacklist')
    requestBlacklist.text = '...'

    subscriptionEmailBlacklist = ET.SubElement(root, 'subscriptionEmailBlacklist')
    subscriptionEmailBlacklist.text = '...'

    user = ET.Comment('<user username="..." password="..." roles="..." />')
    root.append(user)
    return root

#
# add_dataset
#     root: Root DOM for the datasets.xml
#     dataset_name: Name of the dataset
#     file_name: Name of the netCDF file to add to the datasets.xml
#     server_nc_directory: netCDF file path on the ERDDAP server
#     create_virtual_datasets: Whether to combine multiple netCDF files for the
#         same dataset into one virtual ERDDAP dataset.
#
# Add a dataset to the datasets.xml
#
def add_dataset(root, dataset_name, file_name, server_nc_directory, create_virtual_datasets):
    #
    # Root elements
    #
    dataset = ET.SubElement(root, 'dataset', attrib={'type': 'EDDTableFromMultidimNcFiles', 'datasetID': dataset_name, 'active': 'true'})
    reloadEveryNMinutes = ET.SubElement(dataset, 'reloadEveryNMinutes')
    reloadEveryNMinutes.text = '10080'
    updateEveryNMillis = ET.SubElement(dataset, 'updateEveryNMillis')
    updateEveryNMillis.text = '10000'
    fileDir = ET.SubElement(dataset, 'fileDir')
    fileDir.text = server_nc_directory
    fileNameRegex = ET.SubElement(dataset, 'fileNameRegex')
    if create_virtual_datasets:
        fileNameRegex.text = dataset_name + ".*\.nc"
    else:
        fileNameRegex.text = dataset_name + "\.nc"
    recursive = ET.SubElement(dataset, 'recursive')
    recursive.text = 'false'
    pathRegex = ET.SubElement(dataset, 'pathRegex')
    pathRegex.text= '.*'
    metadataFrom = ET.SubElement(dataset, 'metadataFrom')
    metadataFrom.text = 'last'
    preExtractRegex = ET.SubElement(dataset, 'preExtractRegex')
    postExtractRegex = ET.SubElement(dataset, 'postExtractRegex')
    extractRegex = ET.SubElement(dataset, 'extractRegex')
    columnNameForExtract = ET.SubElement(dataset, 'columnNameForExtract')
    removeMVRows = ET.SubElement(dataset, 'removeMVRows')
    removeMVRows.text = 'true'
    sortFilesBySourceNames = ET.SubElement(dataset, 'sortFilesBySourceNames')
    fileTableInMemory = ET.SubElement(dataset, 'fileTableInMemory')
    fileTableInMemory.text = 'false'
    accessibleViaFiles = ET.SubElement(dataset, 'accessibleViaFiles')
    accessibleViaFiles.text = 'false'

    #
    # Get global attributes from netCDF
    #
    nc_file = netCDF4.Dataset(file_name, 'r')
    nc_attrs = nc_file.ncattrs()
    #
    # Top level attributes
    #
    addAttributes = ET.SubElement(dataset, 'addAttributes')
    cdm_data_type = ET.SubElement(addAttributes, 'att', attrib={'name': 'cdm_data_type'})
    cdm_data_type.text = 'Other'
    Conventions = ET.SubElement(addAttributes, 'att', attrib={'name': 'Conventions'})
    Conventions.text = 'COARDS, CF-1.6, ACDD-1.3'
    infoUrl = ET.SubElement(addAttributes, 'att', attrib={'name': 'infoUrl'})
    if 'infourl' in nc_attrs:
        infoUrl.text = nc_file.getncattr('infourl')
    else:
        infoUrl.text = "Local source"
    institution = ET.SubElement(addAttributes, 'att', attrib={'name': 'institution'})
    institution.text = 'US Geological Survey'
    keywords = ET.SubElement(addAttributes, 'att', attrib={'name': 'keywords'})
    # Set the keywords text once the variable names are known
    dataset_license = ET.SubElement(addAttributes, 'att', attrib={'name': 'license'})
    dataset_license.text = 'This USGS product is considered to be in the U.S. public domain'
    sourceUrl = ET.SubElement(addAttributes, 'att', attrib={'name': 'sourceUrl'})
    if 'infourl' in nc_attrs:
        sourceUrl.text = nc_file.getncattr('infourl')
    else:
        sourceUrl.text = 'Local source'
    standard_name_vocabulary = ET.SubElement(addAttributes, 'att', attrib={'name': 'standard_name_vocabulary'})
    standard_name_vocabulary.text = 'CF Standard Name Table v29'
    #subsetVariables = ET.SubElement(addAttributes, 'att', attrib={'name': 'subsetVariables'})
    # Set the keywords text once the variable names are known
    summary = ET.SubElement(addAttributes, 'att', attrib={'name': 'summary'})
    if 'summary' in nc_attrs:
        summary.text = nc_file.getncattr('summary')
    else:
        summary.text = 'No summary info available'
    title = ET.SubElement(addAttributes, 'att', attrib={'name': 'title'})
    if 'title' in nc_attrs:
        title.text = nc_file.getncattr('title')
    else:
        title.text = 'Data from a local source'

    #
    # Variables
    #
    var_names = []
    altitude_and_depth = 'altitude' in nc_file.variables and 'depth' in nc_file.variables

    for nc_var in nc_file.variables:
        var_name = str(nc_var)
        dataVariable = ET.SubElement(dataset, 'dataVariable')
        sourceName = ET.SubElement(dataVariable, 'sourceName')
        sourceName.text = var_name
        destinationName = ET.SubElement(dataVariable, 'destinationName')
        # Special case: ERDDAP is strict about variables named "time""
        # Special case: ERDDAP does not allow a dataset to contain both "altitude" and "depth"
        if (var_name == 'time') or (altitude_and_depth and var_name in ['altitude', 'depth']):
            var_name = 'source_' + var_name
        destinationName.text = var_name
        dataType = ET.SubElement(dataVariable, 'dataType')
        dataType.text = erddap_datatype(nc_file.variables[nc_var].datatype)
        add_field_attrs(dataVariable, nc_var)
        var_names.append(var_name)

    # Now we can set keywords and subsetVariables
    keywords.text = ','.join(var_names)
    #subsetVariables.text = ",".join(var_names)
    
    nc_file.close()

#
# add_field_attrs
#     dataVariable: dataVariable element in the datasets.xml DOM to which to add attributes
#     nc_var: netCDF variable for which to add attributes
#
# Add the appropriate attributes for the given variable.
#
def add_field_attrs(dataVariable, nc_var):
    addAttributes = ET.SubElement(dataVariable, 'addAttributes')

    var_name = str(nc_var)
    ioos_category = get_ioos_category(var_name)
    units = None
    standard_name = None
    if ioos_category == 'Location':
        if var_name.lower() in ['latitude', 'longitude']:
            coordinate_axis_type = var_name[:3].title()
            if coordinate_axis_type =='Lat':
                axis = 'Y'
                units = 'degrees_north'
            else:
                axis = 'X'
                units = 'degrees_east'
        elif var_name.lower() in ['altitude', 'depth']:
            units = 'm'
        standard_name = var_name.lower()

    ioos_category_att = ET.SubElement(addAttributes, 'att', attrib={'name': 'ioos_category'})
    ioos_category_att.text = ioos_category
    long_name_att = ET.SubElement(addAttributes, 'att', attrib={'name': 'long_name'})
    long_name_att.text = var_name.replace('_', ' ').title()
    if standard_name:
        standard_name_att = ET.SubElement(addAttributes, 'att', attrib={'name': 'standard_name'})
        standard_name_att.text = standard_name
    if units:
        units_att = ET.SubElement(addAttributes, 'att', attrib={'name': 'units'})
        units_att.text = units

#
# get_ioos_category
#     var_name: Name of the variable
#
# Determine ioos_category for the given variable.
#
# The current valid values in ERDDAP are Bathymetry, Biology, Bottom Character, Colored Dissolved Organic Matter,
# Contaminants, Currents, Dissolved Nutrients, Dissolved O2, Ecology, Fish Abundance, Fish Species, Heat Flux,
# Hydrology, Ice Distribution, Identifier, Location, Meteorology, Ocean Color, Optical Properties, Other, Pathogens,
# pCO2, Phytoplankton Species, Pressure, Productivity, Quality, Salinity, Sea Level, Statistics, Stream Flow,
# Surface Waves, Taxonomy, Temperature, Time, Total Suspended Matter, Unknown, Wind, Zooplankton Species,
# and Zooplankton Abundance.
#
def get_ioos_category(var_name):
    ret = 'Unknown'
    if var_name.lower() in ['altitude', 'depth', 'latitude', 'longitude']:
        ret = 'Location'
    return ret

#
# erddap_datatype
#     nc_dtype: dtype of netCDF variable
# Returns the ERDDAP mapping of the given netCDF variable
#
def erddap_datatype(nc_dtype):
    ERDDAP_CHAR = 'char'
    ERDDAP_BYTE = 'byte'
    ERDDAP_UBYTE = 'byte'
    ERDDAP_SHORT = 'short'
    ERDDAP_USHORT = 'short'
    ERDDAP_INT = 'int'
    ERDDAP_UINT = 'int'
    ERDDAP_INT64 = 'long'
    ERDDAP_UINT64 = 'long'
    ERDDAP_FLOAT = 'float'
    ERDDAP_DOUBLE = 'double'
    ERDDAP_STRING = 'String'
    ret = None
    if nc_dtype == 'S1':
        ret = ERDDAP_STRING
    elif nc_dtype == 'c':
        ret = ERDDAP_CHAR
    elif nc_dtype == 'i1' or nc_dtype == 'b' or nc_dtype == 'B':
        ret = ERDDAP_BYTE
    elif nc_dtype == 'u1':
        ret = ERDDAP_UBYTE
    elif nc_dtype == 'i2' or nc_dtype == 'h':
        ret = ERDDAP_SHORT
    elif nc_dtype == 'u2':
        ret = ERDDAP_USHORT
    elif nc_dtype == 'i4' or nc_dtype == 'i' or nc_dtype == 'l':
        ret = ERDDAP_INT
    elif nc_dtype == 'u4':
        ret = ERDDAP_UINT
    elif nc_dtype == 'i8':
        ret = ERDDAP_INT64
    elif nc_dtype == 'u8':
        ret = ERDDAP_UINT64
    elif nc_dtype == 'f4' or nc_dtype == 'f':
        ret = ERDDAP_FLOAT
    elif nc_dtype == 'f8' or nc_dtype == 'd':
        ret = ERDDAP_DOUBLE
    return ret


# Extract metadata from returned ScienceBase item records,  
# customize to specs
# could enhance item retrieval from ScienceBase (sb_session)
#
# inputs: filename and partial dictionary from Sciencebase
# outputs: revised metadata dictionary
#
# T. Wellman, BCB Group, USGS
#
def meta_proc(filename, metadata = {}):
 
    attrs = collections.OrderedDict()

    # Collect ScienceBase timestamps ISO formatted UTC
    # ------------------------------------------------
    date_record = []

    logger.debug("\t{}".format('Processed metadata'))

    if 'files' in metadata:
        FileStartDate = ''
        for dfile in metadata['files']:
            if 'dateUploaded' in dfile:
                try:
                    fdate = datetime.strptime(dfile['dateUploaded'],'%Y-%m-%dT%H:%M:%SZ')
                except:
                    try: 
                        fdate = datetime.strptime(dfile['dateUploaded'],'%Y-%m-%dT%H:%M:%S.%fZ')
                    except:
                        fdate = '* timestamp error! *' 
            if fdate != '* timestamp error! *':
                fdate = fdate.strftime('%Y-%m-%dT%H:%M:%SZ')  
                date_record.append(fdate)
            if 'name' in dfile:
                if filename == dfile['name']:
                    FileStartDate = fdate
                elif '.zip' in dfile['name']:
                    if filename.split('.csv')[0] == dfile['name'].split('.zip')[0]:
                        FileStartDate = fdate
        if FileStartDate == '':
            FileStartDate = 'timestamp unavailable *'

    if date_record:
        ItemStartDate = min(date_record)  
    else:
        ItemStartDate = 'timestamp unavailable *'

    # Record current time as ISO formatted UTC
    # -----------------------------------------
    DateNowUtc = datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ') 


    # retrieve + generate global metadata 
    # -------------------------------------
    fmt = '{}\n{}'
    attrs['content'] = fmt.format('OBIS-USA marine biological data consisting of recorded observations', 
                        'of identifiable marine species at a known time and place.')
    if 'title' in metadata:
        attrs['title'] = metadata['title']
    else:
        attrs['title'] = 'Title not yet available'
    if 'body' in metadata:
        fmt = r'(?:\S{81,}|.{1,80})(?:\s+|$)'
        entry = '\n'.join(line.strip() for line in re.findall(fmt, metadata['body']))
        entry = re.sub('<[^>]*>|"', '', entry)
        attrs['summary'] = entry
    else:
        attrs['summary'] = 'Summary information not yet available'
    fmt = '{}\n{}\n{}\n{}'
    attrs['reference'] = fmt.format('The source data is described and accessible in ScienceBase.',                               
           'ScienceBase is a collaborative platform for information sharing and data management',
           'maintained by the U.S. Geological Survey (USGS). Please view item information directly',
           'for additional details: ' + metadata['link']['url'] + ' ') 
    if 'citation' in metadata:
        attrs['citation'] = metadata['citation']
    else:
        attrs['citation'] = 'Journal/report publication not yet available'
    fmt = '{}: {}\n{}: {}\n{}: {}'
    attrs['history'] = fmt.format(ItemStartDate,
           'Phase I - Initial upload of data file(s) to ScienceBase item.', 
           FileStartDate,                               
           'Phase II - Initial QA/QC reprocessing to vocabulary Standards.',
           DateNowUtc, 
           'Phase III - Conversion to netCDF, report diagnostics.')
    agency = 'Biogeographic Characterization Branch, \nCore Science Analytics,'              'Synthesis and Libraries, \nCore Science Systems, \nU.S. Geological Survey' 
    attrs['processor'] = agency

    return attrs


# write dataframe metadata into Xarray dataset
# option - infuse global metadata
# option - infuse variable specific metadata
# option - compare variables to vocabulary standard (dictionary)
#
# T. Wellman, BCB Group, USGS
# --------------------------------------------------------------

def dataframe_metadata(ds, global_metadata, vocab_standard, commands, 
                       err_msg, df_labels, nc_dtypes, vocab_name):                 

    # global metadata
    
    if global_metadata:
        ds.attrs = global_metadata
        
    # variable specfic metadata
    
    var_nconform = 0
    for label in df_labels:
        var_meta = collections.OrderedDict()
        if label in commands['LLAT_specs']:
            var_meta['long_name'] = label.title()
            for key, value in commands['LLAT_specs'][label].items():
                if key != 'destinationName':
                    var_meta[key] = value
                if 'time' in value.lower():
                    if commands['date_convert'].lower() == 'integer':
                        var_meta['units'] = 'days since 1970-01-01'
                    elif commands['date_convert'].lower() == 'datetime':
                        var_meta['units'] = 'datetime64[ns]'
                    else:
                        var_meta['units'] = 'datetime (unknown fmt)'               
        else:
            var_meta['long_name'] = label
            if 'dtype' not in nc_dtypes[label]:
                var_meta['units'] = 'undefined'
            elif 'S' not in nc_dtypes[label]['dtype'] and 'U' not in nc_dtypes[label]['dtype']: 
                    for key,value in commands['variable_units'].items():
                        if key.lower() in label.lower():
                            var_meta['units'] = value
                    if '_date' in label:
                        var_meta['units'] = 'datetime (unknown fmt)'    
                    else:
                        var_meta['units'] = 'unknown' 
            else:
                var_meta['_Encoding'] = commands['string_fmt']
                if any(d in label.lower() for d in ['_date', '_time']):
                    var_meta['units'] = 'datetime (unknown)' 
                
        for key, value in ds[label].attrs.items():
                var_meta[key] = value
        if vocab_standard:
            if label in vocab_standard:
                for key,value in vocab_standard[label].items():
                    entry = '\n'.join(line.strip() for line in re.findall(r'(?:\S{81,}|.{1,80})(?:\s+|$)', value))
                    entry = re.sub('''["']''', '', entry)
                    var_meta[key.lower()] = entry 
            else:
                var_nconform +=1
                var_meta['comment'] = 'non-conforming to ' +  vocab_name
                msg = '\n{}{}{}{}'.format('WARNING: - variable: ',label,', non-conforming to ', vocab_name)
                err_msg.append(msg)

        ds[label].attrs = var_meta


# In[5]:

# brute force (ugly) examination of heterogeneous data types
# examines each dataframe iterator (chunk)
# resolves conflicts between interpreted data chunks
# updates list of dataframe data type (dtype_list) 
#
# T. Wellman, BCB Group, USGS
# -------------------------------------------------------

def data_types(dtype_list, dlist, df, df_labels, commands, err_msg = []): 

    if dtype_list == []:
        for j, d in enumerate(dlist):

            dtype_list.append([])
            dtype_list[j].append(d.base)
            datatype = re.split("([a-zA-Z]+)([0-9]+)", d.name)

            if len(datatype) == 4:
                dtype_list[j].append(datatype[1])
                dtype_list[j].append(int(datatype[2]))
            else:
                dtype_list[j].append(datatype[0])
                dtype_list[j].append('')

            if dtype_list[j][1] == 'object':
                char_max = df[df_labels[j]].dropna().map(len).max()
                dtype_list[j].append(char_max)
            elif dtype_list[j][1] == 'int':
                if commands['netcdf_type'] == 'NETCDF3_CLASSIC':
                    if df[df_labels[j]].dropna().max() > 2147483647:
                        df[df_labels[j]] = df[df_labels[j]].astype('str')
                        char_max = df[df_labels[j]].dropna().map(len).max()
                        dtype_list[j] = [numpy.dtype('O'), 'object', '', char_max]
    else:
        dtype_list2 = []
        for j, d in enumerate(dlist):

            dtype_list2.append([])
            dtype_list2[j].append(d.base)
            datatype = re.split("([a-zA-Z]+)([0-9]+)", d.name)

            if len(datatype) == 4:
                dtype_list2[j].append(datatype[1])
                dtype_list2[j].append(int(datatype[2]))
            else:
                dtype_list2[j].append(datatype[0])
                dtype_list2[j].append('') 

            if dtype_list2[j][1] == 'object':
                char_max = df[df_labels[j]].dropna().map(len).max()
                dtype_list2[j].append(char_max) 
            elif dtype_list2[j][1] == 'int':
                if df[df_labels[j]].dropna().max() > 2147483647:
                    df[df_labels[j]] = df[df_labels[j]].astype('str')
                    char_max = df[df_labels[j]].dropna().map(len).max()
                    dtype_list2[j] = [numpy.dtype('O'), 'object', '', char_max]              

        for j, dtype1 in enumerate(dtype_list):
            if dtype1[0:3] != dtype_list2[j][0:3]:
                fmt = "\t\t{}\n\t\t\t{} {} ({} {})"
                msg = '\n{}{} {} {}'.format('WARNING - column dtypes disagree in chunks: ', 
                                              df_labels[j], dtype1[0], dtype_list2[j][0])
                err_msg.append(msg)
                logger.debug(fmt.format('WARNING - dtypes disagree in chunks: ',
                                 'column label (dtypes): ', df_labels[j], dtype1[0],
                                 dtype_list2[j][0]))

                err = [dtype1[1], dtype_list2[j][1]]
                if err[0] == err[1]:
                    if err[1][2] > err[0][2]:
                        if commands['netcdf_type'] != 'NETCDF3_CLASSIC':
                            dtype_list[j] =  dtype_list2[j]
                            logger.debug("\t\t\t{}".format('action: prescribe larger byte size'))     
                        else:
                            if err[0] != 'int' or err[1][2] <= 32:
                                dtype_list[j] =  dtype_list2[j]
                                logger.debug("\t\t\t{}".format('action: prescribe larger byte size'))       
                elif 'int' in err and 'float' in err:
                    indx = err.index('float')
                    if indx == 1:
                        dtype_list[j] =  dtype_list2[j]
                    logger.debug("\t\t\t{}".format('action: float + integers -> floats'))
                elif 'object' in err:
                    indx = err.index('object')
                    if indx == 1:
                        dtype_list[j] =  dtype_list2[j]
                    logger.debug("\t\t\t{}".format('action: objects + non-objects -> objects'))
            elif dtype_list[j][1] == 'object' and dtype_list2[j][1] == 'object':
                if dtype_list[j][3] < dtype_list2[j][3]:
                    dtype_list[j][3] = dtype_list2[j][3]
                    
    return dtype_list


# In[6]:

# Purpose: csv conversion to NetCDF using (Xarray, Dask) method
#
#   file_name: CSV 
#   source_data_dir: Directory containing CSV source files
#   metadata: adds ScienceBase item metadata + vocabulary metadata
#   erddap_data_dir: Directory in which to write the output netCDF file
#
#  T. Wellman, BCB Group, USGS
#
# approach - lazily processes source csv file using regex to determine formatting,
#            converts specified source CSV to NetCDF with a single "index" (row) dimension,
#            uses xarray-dask dataframe technology for NetCDF conversion + metadata infusion,
#            controls memory consumption using data chunking, incorporates error checking and
#            messaging diagnostics, evaluates dtypes, selects dtype using low-level conflict resolution,
#            rejects source files with critical errrors.
#
# general notes: assumes reasonably well formatted csv, if relevant attempts reformatting from pipe 
# delimited file to comma delimited file, able to pre-filter data using regex procedures constrained
# by prefilter options given in "commands" dictionary.
#
#


def csv_to_nc_dframe(filename, metadata, fpaths, commands, 
                     vocab = [],  p_tasks = [], err_msg = []):
    
    logger.debug("\t{}".format('Implemented dataframe (Xarray, Dask) file conversion method'))
    
    # Regex filter set-up
    # -------------------
    
    regex_filter = re.compile("(%s)" % "|".join(map(re.escape, commands['convert_chars'].keys()))) 
    func = lambda d: commands['convert_chars'][d.string[d.start():d.end()]]

    # Lazily evaluate csv source file, optionally prefilter chars 
    # -----------------------------------------------------------
    
    if not os.path.isfile(fpaths['ncpath']) or commands['proc_overwrite']: 

        with ExitStack() as stack:
            sf = stack.enter_context(open(fpaths['sfpath'], 'r',
                                          encoding= commands['string_fmt'], errors='replace'))

            # read column labels, autodetect delimiters & qouting style,
            # conservatively, process entire file to avoid unknown issues
            # ---------------------------------------------------------
            
            hdr = next(sf)
            delims = [',' , '|' , ';', '\t']
            dcnts = [hdr.count(d) for d in delims]
            max_delim = max(enumerate(dcnts),key=lambda x: x[1])[0]
            delim = delims[max_delim]
            commands['delimiter'] = delim
            splt = re.sub('''["']''', '', hdr).split(delim)
            labels = [h.strip() for h in splt]
            columns = len(labels)
            col_indx =  list(numpy.arange(0, columns))
            
            
            # regex parse format options
            # ----------------------------------
            
            regx_double =  """[\{}]{}""".format(delim,'(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)')
            regx_single =  '''[\{}]{}'''.format(delim,"(?=(?:[^\']*\'[^\']*\')*[^\']*$)")

            # optional prefilter, + auto-activated to force comma delimeter 
            # -------------------------------------------------------------
            
            if commands['prefilter'] or commands['delimiter'] != ',': 
                logger.debug('\t{}'.format('Prefilter active - overwriting source file with prefiltered file'))
                pf = stack.enter_context(open(fpaths['pre_sfpath'],  'a'))
                if hdr.count('"') - hdr.count("'") >= 0:
                    regex_parse = re.compile(regx_double)
                    regx_delm_repl = r'\1,\2'
                else:
                    regex_parse = re.compile(regx_single)
                    regx_delm_repl = r'"\1,\2"'
                regx_delm_srch = r'\s*(\w+),\s*(\w+)'
                hdr = re.sub(regx_delm_srch, regx_delm_repl, hdr)
                read_list = filter_convert_row(hdr, func, regex_filter, regex_parse, col_indx)
                pf.write(','.join(read_list) + '\n')

            
            bulk_stat = collections.defaultdict(int, {'majority' : 0, 'columns': columns, 'rows': 0})     
            qcnt = 0
            
            for indx, read_line in enumerate(sf):

                bulk_stat['rows'] +=1
                qdiff = (read_line.count('"') - read_line.count("'"))
                bulk_stat['majority'] += qdiff

                # optional prefilter, + auto-activated to force comma delimeter 
                # -------------------------------------------------------------
                
                if commands['prefilter'] or commands['delimiter'] != ',':
                    if bulk_stat['majority'] >= 0:
                        regex_parse = re.compile(regx_double)
                    else:
                        regex_parse = re.compile(regx_single)
                    read_line = re.sub(regx_delm_srch, regx_delm_repl, read_line)
                    read_list = filter_convert_row(read_line, func, regex_filter, regex_parse, 
                        col_indx)
                    pf.write(','.join(read_list) + '\n')
                    
                # parse content per line, check for errors
                # -----------------------------------------  
                
                if bulk_stat['majority'] >= 0:
                    sep_line = re.split(regx_double, read_line)
                else:
                    sep_line = re.split(regx_single, read_line)
                if len(sep_line) != columns:
                    logger.critical('\t{}{}'.format('CRITICAL ERROR: - Source file read in row: ', indx))
                    msg = '\n\t{}{}'.format('CRITICAL ERROR: - Source file read in row: ', indx)
                    err_msg.append(msg)
                    commands['QC_pass'] = False
                    return

            # quote formatters
            # -----------------
            
            if bulk_stat['majority'] >= 0:
                q_fmt = ['"', '"{}"', "'"]
                commands['q_fmt'] = q_fmt
                commands['qoute_format'] = '''[\,](?=(?:[^"]*"[^"]*")*[^"]*$)'''
            else:
                q_fmt = ["'", "'{}'", '"']
                commands['q_fmt'] = q_fmt
                commands['qoute_format'] = """[\,](?=(?:[^']*'[^']*')*[^']*$)"""

            if commands['prefilter'] or commands['delimiter'] != ',':
                os.remove(fpaths['sfpath'] )
                os.rename(fpaths['pre_sfpath'], fpaths['sfpath'] )
                commands['delimiter'] = ','
        
            fnt = 0 
            paths = []
            tmpfiles = glob.glob(fpaths['tempdir'] + "*.nc")
            [os.remove(t) for t in tmpfiles] 
            sf.close
            
        # calculate data chunk in terms of table rows 
        # -------------------------------------------
        
        if commands['chunk_elements']:
            commands['chunksize'] =  int(commands['chunk_elements']/columns)
        else:
            commands['chunksize'] = None
        commands['chunksize'] = max(1, commands['chunksize'])
        
        
        # first pass - examine data type consistency
        # reason: dtypes may deviate between data chunks
        # choose data type logically based on chunks
        # ------------------------------------------- 
        
        with open(fpaths['sfpath'], encoding = commands['string_fmt'],
                  errors='replace') as fname:

            dtype_list = []
            dframe = pd.read_csv(fname, chunksize=commands['chunksize'], iterator=True, 
                                 parse_dates=False, quotechar=q_fmt[0],float_precision='round_trip', 
                                 low_memory = False, keep_default_na=False,na_filter=True)
            for df in dframe:

                # remove special characters from variable names 
                # ---------------------------------------------
            
                df.columns = df.columns.str.strip()
                df.columns = df.columns.str.replace(r"[/]", '_per_')
                df.columns = df.columns.str.replace(r"[^A-Za-z0-9]+", '_')                      
                df_labels = df.columns.tolist()
                dlist = df.dtypes.tolist() 
                
                # (re)evaluate data types in datframe iteration
                # ------------------------------------------
                dtype_list = data_types(dtype_list, dlist, df, df_labels, commands, err_msg)
              
            
            # identify time and position variables from LLAT specs
            # and quantititative variables based on name recognititon (terms)
            # -----------------------------------------------------
            
            timevar = [variable for variable in commands['LLAT_specs'] 
                       for key, value in commands['LLAT_specs'][variable].items() if value == 'time']
            posvar = [variable for variable in commands['LLAT_specs'] 
                      for key, value in commands['LLAT_specs'][variable].items() if value != 'time']
            quantvar = [] 
            stringvar = [] 
            for label in df_labels:
                if any(s in label for s in commands['numeric_terms']):
                        quantvar.append(label)
                elif any(s in label for s in commands['string_terms']):
                        stringvar.append(label)
                        

            # assemble NetCDF data type (dtype) dictionary 
            # record largest dimmensions per object
            # include fill value, encoding, and compression
            # ensure LLAT variables, excl. time, are numeric
            # ---------------------------------------------
            
            nc_dtypes = {}
            dtype_dict = {}
            nobj = 0
            for j, datatype in enumerate(dtype_list):
                dtype_dict[df_labels[j]] = datatype[0]
                if datatype[1] == 'object':     
                    if df_labels[j] not in commands['LLAT_specs']:
                        if df_labels[j] not in quantvar:
                            nobj +=1
                            nc_dtypes[df_labels[j]] = {'dtype' : '|S' + str(datatype[3]+2)}      
                            if commands['absent_string'] == '_FillValue':
                                nc_dtypes[df_labels[j]]['_FillValue'] = b''
                # kluge - force data type to string
                elif df_labels[j] in stringvar:
                    nobj +=1
                    nc_dtypes[df_labels[j]] = {'dtype' : '|S32'} 
                if df_labels[j] not in nc_dtypes:
                    nc_dtypes[df_labels[j]] = {} 
                    
                # compress netcdf 
                nc_dtypes[df_labels[j]]['zlib'] = True  
                nc_dtypes[df_labels[j]]['complevel'] = 4
                
            # determine source file qoute style
            # ----------------------------------
            
            if commands['quote_style'] is None:
                elem_total = numpy.multiply(columns, indx+1)
                elem_char_all = numpy.multiply(nobj, indx+1)
                elem_char_obs = numpy.divide(bulk_stat['majority'], 2.0)
                quoted_frac = float(elem_char_obs)/float(elem_char_all)
                qouted_all = float(elem_char_obs)/float(elem_total)
                if qouted_all > 0.990:
                    commands['quote_style'] = csv.QUOTE_ALL
                if quoted_frac > 0.990:
                    commands['quote_style'] = csv.QUOTE_NONNUMERIC
                else:
                    commands['quote_style'] = csv.QUOTE_MINIMAL
   
            paths = []
            fname.seek(0) 
            dframe = pd.read_csv(fname, chunksize=commands['chunksize'], iterator=True,
                                 parse_dates = False, quotechar=q_fmt[0],  
                                 na_filter=True, keep_default_na=False,
                                 float_precision='round_trip', low_memory = False, 
                                 quoting = commands['quote_style'], dtype = dtype_dict,
                                 encoding=commands['string_fmt']) 
            
            for df in dframe:
                  
                logger.debug("\t\t{}{}".format('processing data chunk (rows,columns): ', df.shape))
                
               # modify format (missing values, dates, encoding) 
               # -------------------------------------------------
                                                                 
                # remove special characters from variable names 
                df.columns = df.columns.str.strip()
                df.columns = df.columns.str.replace(r"[/]", '_per_')
                df.columns = df.columns.str.replace(r"[^A-Za-z0-9]+", '_')
                
                if columns == bulk_stat['columns']:
                    for i, col in enumerate(df.columns):
                        if df[col].dtype in ['O', 'S']: 
                            
                            # option, modify fill values
                            if commands['absent_string'] == '_FillValue':
                                df[col] = df[col].fillna(b'') 
                            elif commands['absent_string']:
                                df[col] = df[col].fillna(commands['absent_string']) 
                                
                            # LLAT variables show unexpected type, attempt coerced conversion   
                            if df_labels[i] in posvar:
                                for key, value in commands['LLAT_specs'][df_labels[i]].items():
                                    if key == 'destinationName' and value != 'time':
                                        msg = '\n{}{}{}'.format('WARNING: - variable: ', df_labels[i], ' read as string, forcing conversion to numeric')
                                        err_msg.append(msg)
                                        try:
                                            sample = df[col].loc[~df[col].isnull()].iloc[0]
                                            ssplt = sample.split('.')
                                            if len(ssplt) != 2:
                                                precision = 0
                                            else:
                                                precision = len(sample.split('.')[-1])
                                            entry_ok = df[col].notnull()
                                            df.loc[entry_ok, col] = pd.to_numeric(df[col][entry_ok], errors='coerce')
                                            df[col] = df[col].round(precision)
                                            nc_dtypes[col]['dtype'] = 'float64'
                                        except:
                                            pass
                       
                           # variable is expected numeric, attempt conversion
                            elif df_labels[i] in quantvar: 
                                msg = '\n{}{}{}'.format('WARNING: - variable: ', df_labels[i], ' read as string, forcing conversion to numeric' )
                                err_msg.append(msg)
                                try:
                                    entry_ok = df[col].notnull()
                                    df.loc[entry_ok, col] = pd.to_numeric(df[col][entry_ok], errors='coerce')
                                    nc_dtypes[col]['dtype'] = 'float64'
                                except:
                                    pass
                            
                            # variable is a datetime record
                            elif col in timevar:
                                
                                # attempt to process time variable as utc datetime
                                if commands['date_convert'].lower() == 'datetime':
                                    try:
                                        df[col] = pd.to_datetime(df[col], utc=True)
                                        df[col] = df[col].dt.strftime(commands['date_fmt'])
                                    except ValueError:
                                        pass
                                
                                # attempt to process time variable as utc string data
                                if commands['date_convert'].lower() == 'string':
                                    try:
                                        df[col] = pd.to_datetime(df[col], utc=True)
                                        df[col] = df[col].dt.strftime(commands['date_fmt'])
                                        nc_dtypes[col] = {'dtype' : 'S22'} 
                                    except ValueError:
                                        pass
                                    
                                # attempt to process time variable as integer-days since (1970) basedate
                                elif commands['date_convert'].lower() == 'integer':
                                    try:
                                        df[col] = pd.to_datetime(df[col], utc=True)
                                        base = pd.to_datetime('1970-01-01', utc=True)
                                        df[col] = (df[col] - base).dt.days 
                                        df[col] = df[col].astype(int)
                                        nc_dtypes[col] = {'dtype' : 'int'}
                                    except ValueError:
                                        pass          
                                    
                            else:
                                # re-encode dtype so specified string length is honored
                                # between chunks ... xarray interpretive issue.
                                df[col] = df[col].str.encode('ascii', errors = 'ignore')
                                df[col] = df[col].str.decode(commands['string_fmt'])
                                
                        # variable should be converted to string 
                        elif df_labels[i] in stringvar: 
                                msg = '\n{}{}{}'.format('WARNING: - variable: ', df_labels[i], 
                                                        ' read as numeric, forcing conversion to string' )
                                err_msg.append(msg)
                                try:
                                    df[col] = df[col].astype(str)
                                except:
                                    pass 
                                
                else:
                    logger.critical('\t{}'.format('CRITICAL ERROR: - disagreement in data columns'))
                    msg = '\n\t{}'.format('CRITICAL ERROR: - disagreement in data columns')
                    err_msg.append(msg)
                    commands['QC_pass'] = False
                    return
               
                # convert dataframe to xarray dataset
                # -------------------------------------

                ds = xr.Dataset.from_dataframe(df)

                # save to netCDF, as single or chunked aggregate process
                # ------------------------------------------------------

                if not commands['chunksize'] or bulk_stat['rows'] <= commands['chunksize']:
                    
                    dataframe_metadata(ds, metadata, vocab, commands, err_msg, df_labels,  
                                               nc_dtypes, commands['vocab_name']) 
                    ds.to_netcdf(fpaths['ncpath'], mode = 'w', 
                                 format = commands['netcdf_type'], encoding = nc_dtypes) 
                    commands['QC_pass'] = True
                    ds.close()
                else:
                    
                    fnt+=1
                    chunkfile = fpaths['tempdir'] + 'netcdf_chunk_' + str(fnt) + '.nc'
                    ds = ds.chunk(chunks=commands['chunksize'], name_prefix='xarray-', 
                                  token=None, lock=True)
                    ds.to_netcdf(chunkfile, mode = 'w', 
                                 format = commands['netcdf_type'], encoding = nc_dtypes)
                    paths.append(chunkfile)
                    ds.close()
                    
            
            # if chunked, attempt to re-assemble via Dask processing through Xarray
            # optionally, attempt to enforce Climate-Forceast compliancy (may fail)
            # ------------------------------------------------------------------
            
            if commands['chunksize']: 
                if bulk_stat['rows'] > commands['chunksize']:
                    
                    commands['QC_pass'] = True
                    decodeCF = True
                
                    if commands['cf_comply']:
                        try:
                            ds = xr.open_mfdataset(paths, decode_cf= decodeCF, concat_dim='index')
                            ds = ds.chunk(chunks=commands['chunksize'], name_prefix='xarray-', 
                                          token=None, lock=True)
                            dataframe_metadata(ds, metadata, vocab, commands, err_msg, 
                                               df_labels, nc_dtypes, commands['vocab_name'])
                            ds.to_netcdf(fpaths['ncpath'], mode = 'w',
                                         format = commands['netcdf_type']) 
                            ds.close()
                        except:
                            decodeCF = False    
                            
                        if not decodeCF:
                            try:
                                ds = xr.open_mfdataset(paths, decode_cf= decodeCF, concat_dim='index') 
                                ds = ds.chunk(chunks=commands['chunksize'], name_prefix='xarray-', 
                                              token=None, lock=True)
                                dataframe_metadata(ds, metadata, vocab, commands, err_msg, 
                                                   df_labels, nc_dtypes, commands['vocab_name'])
                                ds.to_netcdf(fpaths['ncpath'], mode = 'w', 
                                             format = commands['netcdf_type']) 
                                ds.close()
                            except:
                                commands['QC_pass'] = False
                    else:
                        
                        decodeCF = False
                        try:
                            ds = xr.open_mfdataset(paths, decode_cf= decodeCF, concat_dim='index') 
                            ds = ds.chunk(chunks=commands['chunksize'], name_prefix='xarray-',
                                          token=None, lock=True)
                            dataframe_metadata(ds, metadata, vocab, commands, err_msg, 
                                               df_labels, nc_dtypes, commands['vocab_name'])
                            ds.to_netcdf(fpaths['ncpath'], mode = 'w', format=commands['netcdf_type'])
                            ds.close()
                        except:
                            commands['QC_pass'] = False
                    
                    if commands['QC_pass']:
                        logger.debug('\t\t{}{}'.format('climate-Forcast compliancy (Xarray) - SET TO: ', decodeCF))
                        msg = '\n{}{}\n'.format('Climate-Forcast compliancy (Xarray) - SET TO: ', decodeCF)
                        err_msg.append(msg)
                    else:
                        logger.debug('\t{}'.format('CRITICAL ERROR: - conversion to netCDF FAILED'))
                        msg = '\n\t{}'.format('CRITICAL ERROR: conversion to netCDF FAILED')
                        err_msg.append(msg)
                        
                        if os.path.isfile(fpaths['ncpath']):
                            os.remove(fpaths['ncpath'])
                        fileQCfail = './QC_FAIL/' + fpaths['sfpath'].split('/')[-1] 
                        os.rename(fpaths['sfpath'], fileQCfail) 
                        [os.remove(p) for p in paths]

            if not commands['QC_pass']:
                errmsg = "** QC/PROCESSING FAIL ** - during conversion, sending files to ./QC_FAIL directory"
                logger.warning("\n\t{}\n".format(errmsg))
                p_tasks.append(errmsg)
                dumpfiles = [fpaths['sfpath'], fpaths['ncpath'],fpaths['rconpath']]
                for orig_file in dumpfiles:
                    if os.path.isfile(orig_file):
                        mvfile = './QC_FAIL/' + orig_file.split('/')[-1] 
                        os.rename(orig_file, mvfile)
            else:
                logger.debug('\t\t{}'.format('... completed file conversion'))
                p_tasks.append('processed source csv')
                p_tasks.append('created netcdf')
                
            return
            
    elif os.path.isfile(fpaths['ncpath']):
        commands['QC_pass'] = True
        return
        


# In[7]:

# Purpose: initialize/modify default input parameters
# accepts either commandline options or *yaml file
#

def arg_overwrite(opts, args, commands):
    
    # optional default overwrite - *.yaml config. file  - T. Wellman
    # -------------------------------------------------

    if any(y in arg for arg in args for y in ['.yml', '.yaml']):
        if 'yaml' in sys.modules:
            with open(args[0], 'r') as options:
                try:
                    opt_dict = yaml.load(options)
                except yaml.YAMLError as exc:
                    print(exc)
                for key, val in opt_dict.items():
                    if key in commands:
                        commands[key] = val
                    else:
                        logger.debug('{} - yaml option was not identified'.format(key))
        else:
            logger.critical('{}\n'.format('Please install yaml package or use default getopt command format'))
            sys.exit()
            

    # optional getopt overwrite - command line inputs
    # -----------------------------------------------   

    elif 'opts' in locals():
        for o, a in opts:
            if o in '--workdir':
                commands['workdir'] = a
            elif o in '--collectionid':
                commands['collection_id'] = a
            elif o in '--sourcedir':
                commands['source_data_dir'] = a
            elif o in '--erddapdir':
                commands['erddap_data_dir'] = a
            elif o in '--serverdir':
                commands['server_nc_directory'] = a
            elif o in '--file_overwrite':
                commands['file_overwrite'] = a
            elif o in '--proc_overwrite':
                commands['proc_overwrite'] = a
            elif o in '--convert_method':
                commands['convert_method'] = a 
            elif o in '--recon_data_dir':
                commands['recon_data_dir'] = a
            elif o in '--report_dir':
                commands['report_dir'] = a
            elif o in '--tempdir':
                commands['tempdir'] = a
            elif o in '--compare_csv2csv':
                commands['compare_csv2csv'] = a
            elif o in '--prefilter':
                commands['prefilter'] = a
            elif o in '--postfilter':
                commands['postfilter'] = a
            elif o in '--dump_csv':
                commands['dump_csv'] = a
            elif o in '--error_report':
                commands['error_report'] = a
            elif o in '--max_report':
                commands['max_report'] = int(a)
            elif o in '--table_output':
                commands['table_output'] = a
            elif o in '--int_float_accept':
                commands['int_float_accept'] = a
            elif o in '--date_fmt':
                commands['date_fmt'] = a
            elif o in '--date_convert':
                commands['date_convert'] = a
            elif o in '--string_fmt':
                commands['string_fmt'] = a
            elif o in '--csv_qoute_logic':
                commands['csv_qoute_logic'] = a
            elif o in '--chunksize':
                commands['chunk_elements'] = int(a)
            elif o in '--cf_comply':
                commands['cf_comply'] = a
            elif o in '--absent_string':
                commands['absent_string'] = a
            elif o in '--netcdf_type':
                commands['netcdf_type'] = a 
            elif o in '--fname_ext':
                commands['fname_ext'] = a
            elif o in '--virtual_datasets':
                commands['create_virtual_datasets'] = True
            elif o in '--window':
                commands['window'] = int(a)
            elif o in '--rowsperfile':
                commands['rows_per_file'] = int(a)
            elif o in '--sample':
                commands['sample_size'] = int(a)
            elif o in '--verbose':
                commands['verbose'] = True
            elif o in '--only_csvs':
                commands['fetch_metadata'] = False
                commands['fetch_csvs'] = True
                commands['create_netcdf_files'] = False
                commands['create_datasets_xml'] = False
            elif o in '--only_netcdf':
                commands['fetch_metadata'] = True
                commands['fetch_csvs'] = False
                commands['create_netcdf_files'] = True
                commands['create_datasets_xml'] = False
            elif o in '--only_datasets_xml':
                commands['fetch_metadata'] = False
                commands['fetch_csvs'] = False
                commands['create_netcdf_files'] = False
                commands['create_datasets_xml'] = True
            elif o in '--help':
                get_ipython().system(u'usage()')
                exit(0)
            else:
                assert False, 'unhandled option'
                
    return commands


# In[8]:

#
# routine to create/clean work directories 
# ------------------------------------------------------
#
# T. Wellman, BCB Group, USGS

def set_directories(commands):
    if commands['workdir']:
        os.chdir(commands['workdir'])
    folder_list = ['source_data_dir', 'erddap_data_dir', 'server_nc_directory',
                   'recon_data_dir', 'report_dir', 'tempdir']
    for folder in folder_list:
        if not os.path.isdir(commands[folder]):
            os.makedirs(commands[folder])
        else:
            commands[folder] 
            srch = (commands[folder] + '/.*').replace('//','/')
            files = glob.glob(srch)
            for f in files:
                os.remove(f)
    if not os.path.isdir('./QC_FAIL'):
            os.makedirs('./QC_FAIL') 
    else:
        files = glob.glob('./QC_FAIL/.*')
        for f in files:
            os.remove(f)
    commands['prefilter_dir'] = './prefilter_data'
    if not os.path.isdir(commands['prefilter_dir']):
        os.makedirs(commands['prefilter_dir'])
    else:
        srch = (commands['prefilter_dir'] + '/.*').replace('//','/')
        files = glob.glob(srch)
        for f in files:
            os.remove(f)


# In[9]:

#
# Routine to activate logger function (default output is to file)
# optional print to console, sets log print level (DEBUG, INFO, WARNING etc.)
#
# T. Wellman, BCB Group, USGS

# create new message logger (100 level)
# -------------------------------------
def add_logger(self, message, *args, **kws):
    self.log(100, message, *args, **kws)


# activate logger with opts
# -------------------------------  
def logform(screen, log_level):
    
    global logger
    
    logging.addLevelName(100, 'MESSAGE')
    logging.Logger.message = add_logger
    
    # report log results to file
    # ----------------------------------
    logging.basicConfig(level=log_level,
                        format= '%(asctime)s.%(msecs)03d %(name)-12s:  %(levelname)-8s %(message)s', 
                        datefmt= '%Y-%m-%d %H:%M:%S',
                        filename='obis_processor.log',
                        filemode='w')         
    logger = logging.getLogger(__name__)
    logger.setLevel(log_level)
    if logger.hasHandlers():
        logger.handlers.clear()
    
    # option, also report logs to console 
    # -----------------------------------
    if screen:
        console = logging.StreamHandler()
        logger.addHandler(console)
        formatter = logging.Formatter(
            fmt = '%(asctime)s.%(msecs)03d - %(levelname)8s - %(message)s', 
            datefmt= '%Y-%m-%d %H:%M:%S')
        console.setFormatter(formatter)
        console.setLevel(log_level)
        
    return


# In[10]:

# setup file paths (per data item)
# remove outdated files, if relevant
# -----------------------------------
#
# T. Wellman, BCB Group, USGS

def setup_fpaths(filename, commands, p_tasks = [], err_msg = []):
        
    fpaths = {}
    fpaths['tempdir']  = commands['tempdir'] 
    fpaths['sfpath']   = os.path.join(commands['source_data_dir'], filename)
    fpaths['ncpath']   = os.path.join(commands['erddap_data_dir'], filename).rsplit('.csv')[0] + '.nc' 
    fpaths['rconpath'] = os.path.join(commands['recon_data_dir'], 
                                           filename).rsplit('.csv')[0] + commands['fname_ext'] + '.csv' 
    fpaths['cmp_report'] = os.path.join(commands['report_dir'], 
                                             filename).rsplit('.csv')[0] + '_compareLOG.txt'
    fpaths['var_report'] = os.path.join(commands['report_dir'], 
                                             filename).rsplit('.csv')[0] + '_variableINFO.txt' 
    
    
    if os.path.exists(fpaths['sfpath']):
        p_tasks.append('Identified source file in local directory')
        
        # eliminate outdated files
        sf_tstamp = os.path.getmtime(fpaths['sfpath'])
        files = [fpaths['ncpath'], fpaths['rconpath'], fpaths['cmp_report'], fpaths['var_report']]
        for file in files:
            if os.path.exists(file):
                pf_tstamp = os.path.getmtime(file)
                if sf_tstamp > pf_tstamp:
                    os.remove(file)
        
    if commands['prefilter']:
        p_tasks.append('Prefiltered source data')
        err_msg.append("\n{}\n{}\n\n".format( 'NOTE: pre-filtering is active (pre-filtered source file)',
                        '   results may be different than using raw source data'))
    prefile = filename.rsplit('.csv')[0] + '_prefltr' + '.csv'
    fpaths['pre_sfpath'] = os.path.join(commands['prefilter_dir'], prefile)
    os.remove(fpaths['pre_sfpath']) if os.path.exists(fpaths['pre_sfpath']) else None
    return fpaths
        


# In[11]:

# Retrieve and process vocabulary standards from
# Biodiversity Information Standards (TDWG) and ESIP SPARQL endpoint.

def retrieve_vocab():

    # parse darwin core terms from html
    # -----------------------------------
    response = request.urlopen('http://rs.tdwg.org/dwc/terms/index.htm') 
    bs = BeautifulSoup(response, "lxml")


    # grab D. Core terminology table in adhoc way
    # TDWG html tables currently have no id/name
    # -------------------------------------------
    bs_table = bs.findAll('table')[4]


    # Process dictionary of Darwin Core terms:
    # -----------------------------------------
    vocab_terms = collections.OrderedDict()
    for tbody in bs_table.findAll('tbody'):
        for tr in tbody.findAll('tr'):

            # process term name
            ta = tr.find('a')
            if ta:
                if 'Term Name' in ta.string:
                    key, term_name = ta.string.split(':')[-2:]
                    term_name = term_name.strip()
                    vocab_terms[term_name] = collections.OrderedDict()

            # process term information + links
            td = tr.findAll('td')
            if td:
                try:
                    label = td[0].text.replace(':','').lower()
                    vocab_terms[term_name][label] = td[1].text + ' '
                except:
                    pass
                href = td[1].find('a')

                # replace details link with url
                if href and td[0].text == 'Details:':
                    dir_base = 'http://rs.tdwg.org/dwc/terms'
                    dir_detail = href.attrs['href'].rsplit('./')[-1]
                    url = os.path.join(dir_base, dir_detail)
                    vocab_terms[term_name]['details'] = url
                    
    # from query results, determine key order (lowercase) and record length
    DCkey_order = list(list(vocab_terms.items())[0:1][0][1].keys())
    DCkey_order = [ dc.lower() for dc in DCkey_order ]
    orig_dic_length = len(vocab_terms.keys())


    # Next, retrieve OBIS vocabulary standards 
    # from ESIP SPARQL endpoint using triple stores
    # ---------------------------------------------

    query_tag = ''' PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    SELECT ?label ?identifier ?section ?definition ?details 
    WHERE {
        ?identifier ?pred1 ?definition; ?pred2 ?label; ?pred3 ?details; ?pred4 ?section  .
        FILTER(regex(str(?identifier), "http://mmisw.org/ont/ioos/marine_biogeography" )) . 
        FILTER(regex(str(?pred1),"/Defintion")) .
        FILTER(regex(str(?pred2),"/Term")) .
        FILTER(regex(str(?pred3),"/Reference")) .
        FILTER(regex(str(?pred4),"/Section")) .
        }'''

    sparql = SPARQLWrapper("http://cor.esipfed.org/sparql")
    sparql.setQuery(query_tag)
    sparql.setReturnFormat(JSON)

    results = {}
    try :
        results = sparql.query()
        results = results.convert()
    except:
        print('exception occurred - reconfigure query')
    try:
        proc_list = results['results']['bindings']
    except:
        proc_list = {}
        print('no vocab data - reconfigure query')


    # Combine unique search terms from both queries 
    # into a single record dictionary 
    # -------------------------------

    for dic in proc_list:
        if 'label' in dic:
            if isinstance(dic['label'],dict):
                label = dic['label']['value'].lower()
            else:
                label = dic['label'].lower()
            if label not in vocab_terms.keys():
                vocab_terms[label] = {}
                otherkeys = [key.lower().lower() for key,value in dic.items() if key.lower() != 'label'] 
                variable_dic = collections.OrderedDict()
                for key in DCkey_order:
                    if key in otherkeys:
                        if isinstance(dic[key],dict):
                            variable_dic[key] = dic[key]['value'].lower()
                        else:
                            variable_dic[key] = dic[key].lower()
                    else:
                        variable_dic[key] = ''
                for key in otherkeys:
                    if key not in variable_dic:
                        if isinstance(dic[key],dict):
                            variable_dic[key] = dic[key]['value'].lower()
                        else:
                            variable_dic[key] = dic[key].lower()
                vocab_terms[label] =   variable_dic     

    # from other query results, determine additional vocab terms from ESIP sparql endpoint
    add_terms = len(vocab_terms.keys()) - orig_dic_length 
    print("\n\t{}{}\n\t{}{}".format('TDWG vocab terms: ', orig_dic_length, 'ESIP vocab terms: ', add_terms))
    
    return DCkey_order, vocab_terms  


# Definitions to Query ESIP-USGS SPARQL Marine Vocabulary terms 
# return as *.json result, returns None if unsuccessful
# incorporates low-level fuzzy queries

def esip_srch_string(opt, vname):
    
    # base sparql ESIP query 
    query_base = """
    SELECT ?label ?identifier ?section ?definition ?details 
    WHERE { ?identifier ?pred1 ?definition; ?pred2 ?label; ?pred3 ?details; ?pred4 ?section  .
        ?pred1 rdfs:label ?predobj
        FILTER(regex(str(?pred1),"/Defintion")) .
        FILTER(regex(str(?pred2),"/Term")) .
        FILTER(regex(str(?pred3),"/Reference")) .
        FILTER(regex(str(?pred4),"/Section")) .
        FILTER(regex(str(?identifier), "ioos" ,"i" )) . 
        """

    # sparql basic fuzzy search logic
    if opt == 'verbatim_case':
        add_filter = """{}"{}"{}""".format( "FILTER(str(?label) = ", vname, '). \n}')  
        query_tag = query_base + add_filter
    if opt == 'verbatim_no_case': 
        add_filter = """{}"{}"{}""".format( "FILTER(lcase(str(?label)) = ", vname.lower(), '). \n}')  
        query_tag = query_base + add_filter
    elif opt == 'no_special_chars': 
        vname_regx = re.sub('[^A-Za-z0-9]+', '', vname)
        add_filter = """{}"{}"{}""".format( "FILTER(lcase(str(?label)) = ", vname_regx.lower(), '). \n}')   
        query_tag = query_base + add_filter
    elif opt == 'split_terms':
        vname_regx = re.sub('[^A-Za-z0-9]+', '_', vname)
        vname_regx = re.split(r'[_]', vname_regx)
        add_filter = ''
        for v in vname_regx:
            add_filter = add_filter + """{}"{}"{}""".format( "FILTER(regex(str(?label), ", v, ', "i" )).\n') 
        add_filter = add_filter + '}' 
        query_tag = query_base + add_filter
    else:
        query_tag = None
        
    return query_tag
    
def sparql_query(query_tag, sparql_endpoint):
    
    sparql = SPARQLWrapper(sparql_endpoint)
    sparql.setQuery(query_tag)
    sparql.setReturnFormat(JSON)
    results = {}
    
    try :
        results = sparql.query()
        results = results.convert()
    except:
        print('No data returned - reconfigure query')
    
    if not results:
        results = None
    elif 'results' in results:
        if 'bindings' in results['results']:
            if not results['results']['bindings']:
                results = None
        else:
            results = None
                
    return results


# Definition to retrieve and process vocabulary standards from
# Biodiversity Information Standards (TDWG) and ESIP SPARQL endpoint.
# Perform match tests using case insensitive matches and fuzzy matches
# (removing special chars and split word searches (split on spec. chars).
# Considers unique search result valid.
# Identifies non-unique results where applicable.
# Attempts update of missing NetCDF variable information.
#
# TDWG vocab match attempts are case sensitive and insensitive (w/o fuzzy)

def eval_vocab(commands):

    #
    # vocab search path and parameters
    #
    ncfiles = glob.glob(commands['erddap_data_dir'] + "/*.nc")
    sparql_endpoint = "http://cor.esipfed.org/sparql" 
    qtypes = ['no_special_chars', 'split_terms']
    msg = 'information from Biodiversity Information Standards (TDWG) and USGS-ESIP SPARQL endpoint'
                               
    #
    # process local NetCDF, determine variable matches where information is missing
    #
    DCkey_order, vocab_terms = retrieve_vocab()
    vocab_terms_lower = [key.lower() for key in vocab_terms.keys()]

    for ncfile in ncfiles:
        
        # proceed when proc_overwrite is active or vocab file is not present in folder
        dcfname = os.path.join(commands['report_dir'], ncfile.split('/')[-1]).rsplit('.nc')[0] + '_variableINFO.json' 
        if commands['proc_overwrite'] or not os.path.isfile(dcfname):
            
            # track missing variable information, match type, and match candidates
            match_types = collections.OrderedDict([('verbatim', 0),
                                                   ('ignore_case', 0),
                                                   ('fuzzy', 0),
                                                   ('nonunique', 0),
                                                   ('no_match', 0),])
            srch_dic = collections.OrderedDict([('explanation', 'information on netcdf variables and match types'),
                                        ('missing_vars', []),
                                        ('match_count' , match_types),
                                        ('possible_matches' , {})])

            # loop through NetCDF variables, attempt to find missing variable information
            
            with netCDF4.Dataset(ncfile, 'r+') as ncc:
        
                dim_vars = [d for d in ncc.dimensions.keys() if 'string' not in d ]

                for key in ncc.variables.keys():
                    if key not in dim_vars:

                        # matches exactly 
                        if key in vocab_terms.keys():
                            srch_dic['match_count']['verbatim'] += 1
                            for sdkey in vocab_terms[key].keys(): 
                                ncc.variables[key].setncattr(sdkey, vocab_terms[key][sdkey]) 
                            if 'comment' not in vocab_terms[key].keys():
                                ncc.variables[key].setncattr('comment', msg)

                        # matches as case insensitive
                        elif key.lower() in vocab_terms_lower:
                            srch_dic['match_count']['ignore_case'] += 1
                            for sdkey in vocab_terms[key.lower()].keys():
                                ncc.variables[key].setncattr(sdkey, vocab_terms[key.lower()][sdkey]) 
                            if 'comment' not in vocab_terms[key.lower()].keys():
                                ncc.variables[key].setncattr('comment', msg)

                        # attempt term association using low-level fuzzy logic
                        else:
                            hits = 0; i = 0; ideal_hits = 0
                            while hits != 1 and i < len(qtypes):
                                query_tag = esip_srch_string(qtypes[i], key)
                                if query_tag is not None:
                                    results = sparql_query(query_tag, sparql_endpoint)
                                    if results is not None:
                                        hits = len(results['results']['bindings'])
                                        if hits == 1 or hits > ideal_hits:
                                            ideal_hits = hits
                                            best_results = results['results']['bindings']       
                                i += 1
                            if ideal_hits == 1:
                                srch_dic['match_count']['fuzzy'] += 1
                                best_results = best_results[0]
                                label = best_results['label']['value']
                                vocab_terms[label] = {}
                                otherkeys = [okey.lower().lower() for okey,value in best_results.items()                                              if okey.lower() != 'label'] 
                                variable_dic = collections.OrderedDict()
                                for varkey in DCkey_order:
                                    if varkey in otherkeys:
                                        variable_dic[varkey] = best_results[varkey]['value']
                                        ncc.variables[key].setncattr(varkey, variable_dic[varkey])
                                if 'label' in best_results:
                                    variable_dic['label'] = label 
                                    ncc.variables[key].setncattr('alias_name', label)
                                for varkey in otherkeys:
                                    if varkey not in variable_dic:
                                        variable_dic[varkey] = best_results[varkey]['value']
                                        ncc.variables[key].setncattr(varkey, variable_dic[varkey])
                                if 'comment' not in variable_dic.keys():
                                    ncc.variables[key].setncattr('comment', msg)
                                    variable_dic['comment'] = msg
                                vocab_terms[label] = variable_dic   
                            elif ideal_hits > 1:
                                srch_dic['match_count']['nonunique'] += 1
                                srch_dic['possible_matches'][key] = best_results
                                ncc.variables[key].setncattr('comment', 'Refer to ScienceBase record for explanation') 
                            elif ideal_hits == 0:
                                srch_dic['match_count']['no_match'] += 1
                                srch_dic['missing_vars'].append(key)
                                ncc.variables[key].setncattr('comment', 'Refer to ScienceBase record for explanation') 

                # output vocabulary match results               
                with open(dcfname, 'w', encoding = 'utf-8', errors = 'replace') as outfile:
                    json.dump(srch_dic, outfile, indent=4)
        
    
    # Create combined missing vocabulary file from individual vocab files
    logger.debug('** Writing file of global missing label definitions (missing_vocab.json)')
       
    varfiles = glob.glob(commands['erddap_data_dir'] + "/*_variableINFO.json")
    missing_vocab = collections.OrderedDict([('explanation', 'bulk accounting of all unassociated netcdf variables'),
                                        ('missing_vars', [])])
    for vfile in varfiles:
         with open(vfile, 'r') as infile:
            jdic = json.load(infile)
            for key in jdic['missing_vars']:
                if key not in missing_vocab:
                    missing_vocab['missing_vars'].append(key)
                    
    with open('missing_vocab.json', 'w', encoding = 'utf-8', errors = 'replace') as outfile:
        json.dump(missing_vocab, outfile, indent=4)
        
    # output revised vocabulary file
    datenow = str(datetime.utcnow().strftime('%Y-%m-%dT%H-%M-%SZ'))
    dcfname = ('NetCDF_revised_vocab' + datenow).replace(' ','_') + '.json'
    with open(dcfname, 'w', encoding = 'utf-8', errors = 'replace') as outfile:
        json.dump(vocab_terms, outfile, indent=4)


# In[12]:


#######################################################################
#                      main processing section
#######################################################################

# processor is used to process list of ScienceBase records and local files,
# creates netCDF files from *.csv files, performs QA/QC, and infuses metadata.
# Uses either messytables or dataframe (xarray, dask) approach for processing data files. 
# Includes process reporting, content filtering, error messages, and ERDDAP datasets.xml file generation.


# if applicable, suppress print-to-screen from imported modules
# -------------------------------------------------------------
def supress_print(func):
    def wrapper(*args, **kargs):
        with open(os.devnull, 'w') as devnull:
            with redirect_stdout(devnull):
                func(*args, **kargs)        
    return wrapper


@supress_print
def run(**kwargs):

    # update commands dictionary, activate logger  
    # ----------------------------------------------
    
    commands = default_inputs()
    for c in kwargs:
        if c in commands:
            commands[c] = kwargs[c]
    logform(commands['log_screen'], commands['log_level'])
    logger.message('** Processor activated **')
            
            
    # create + clean local directories
    # --------------------------------
    set_directories(commands)
        
            
    # optionally, retrieve files and source metadata (locally, ScienceBase) 
    # ------------------------------------------------------------------
    sb = SbSession()
    listdir, sources = retrieve_source_files(sb, commands)
    if commands['proc_limit']:
        listdir = listdir[0:min(len(listdir), commands['proc_limit'])]
    
                                
    # optionally retrieve specified vocabulary definitions
    # -------------------------------------------------------
    if commands['vocab']:
        dc = open(commands['vocab'], 'r', 
                  encoding= commands['string_fmt'], errors='replace')
        vocab = json.load(dc)
    else:
        vocab = []
        
        
    # process file-by-file 
    # -----------------------
    convert_chars =  {k.lower(): v for k, v in commands['convert_chars'].items()}
    logger.debug("{}".format('** Processing available data files'))
    for k, filename in enumerate(listdir):
        
        # initialize messaging
        # ---------------------
        err_msg = []
        p_tasks = []
        logger.debug("{}{}{}{}".format('Processing file # ', k+1,': ', filename))
        msg = '\n\nProcessing commands: \n'
        for key, value in commands.items():
            msg = "{} key: {} value: {} \n".format(msg, key, value) 
        err_msg.append(msg + '\n\n')
        
        # setup file paths, purge if updating source file
        # ------------------------------------------------
        fpaths = setup_fpaths(filename, commands, p_tasks, err_msg)
            
        
        # 1) optionally, generate *.netCDF file from source *.csv file
        # ------------------------------------------------------------
        if commands['create_netcdf_files']:
            
            # process ScienceBase item metadata
            if filename in sources:
                item_metadata = sources[filename]
                cmp_metadata = meta_proc(filename, item_metadata)
            else:
                cmp_metadata = {}
                

            # file conversion (uses "dataframe" or "messytables" approach)
            if commands['convert_method'] == "dataframe":
                csv_to_nc_dframe(filename, cmp_metadata, fpaths, commands,
                                 vocab, p_tasks, err_msg)
            else: 
                csv_to_nc_messy(filename, commands['source_data_dir'], 
                                cmp_metadata, commands) 
  
        
        # 2) optionally evaluate file conversion, w/ filtering + QC reporting 
        # -------------------------------------------------------------------
        if commands['compare_csv2csv']:
            
            if commands['QC_pass']: 
                
                # reconstruct csv files from netCDF
                regenerate_csv(fpaths, commands, err_msg, p_tasks) 

                # compare original (source) csv to regenerated csv
                compare_csv(fpaths, commands, err_msg, p_tasks)
                
        
        # 3) output qc/processing tasks (other reporting in step 2):
        # -----------------------------------------------------------
        msg = str(("\t{}".format('Tasks completed:')))
        for p in p_tasks:
            logger.debug("\t{}".format(p))     

            
    # 4) optionally, evaluate NetCDF variable descriptions, find missing variable metdata, 
    #    use online vocabulary resources to perform verbatim and fuzzy searches
    # ----------------------------------------------------------------------
    if commands['vocab']:
         logger.debug('{}'.format("NetCDF metadata check skipped - use 'python -m obis.metadata_check' via command line"))        
            
            
    # 5) optionally create the ERDDAP datasets.xml file using netCDF files
    # ---------------------------------------------------------------------
    if commands['create_datasets_xml'] :
        write_datasets_xml(commands, "./")
        logger.debug('   {}'.format('datasets.xml file created'))
    
    logger.message("{}\n\n".format('** Processor terminated **'))

    
# **** main program ****


logform(True, logging.INFO) 

if __name__ == "__main__":
    
    opt_fields = [
    'collectionid=', 'sourcedir=', 'erddapdir=', 'serverdir=', 'only_csv', 'only_netcdf',
    'only_datasets_xml', 'virtual_datasets', 'window=', 'rowsperfile=', 'file_overwrite', 'proc_overwrite',
    'sample=', 'verbose', 'help', 'recon_data_dir=', 'report_dir=', 'tempdir=', 'compare_csv2csv',
    'error_report', 'postfilter', 'dump_csv', 'max_report=', 'table_output', 'convert_method=',
    'chunksize=', 'prefilter', 'date_convert', 'date_fmt=', 'int_float_accept']

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'fc', opt_fields)
    except getopt.GetoptError as err:
        logger.error(str(err))
        get_ipython().system(u'usage(); exit(2) ')
    
    # retrieve processing commands
    # -----------------------------
    commands = default_inputs()    
    commands = arg_overwrite(opts, args, commands)
    
    run(**commands)
                       


# In[13]:

# !ncdump -h '/Users/twellman/erddap_data/nc_test/DryTortugasReefVisualCensus2004_Event.nc'
# !ncdump -h '/Users/twellman/erddap_data/nc_test/DryTortugasReefVisualCensus2004_measurementOrFact.nc'


# In[ ]:



