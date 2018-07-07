
# coding: utf-8

# In[1]:


# Script to modify metadata variable attribute records
# from a directory of NetCDF.
#
# Performs vocabulary retrieval and simple processing
# Extracts and processes vocabulary data from 
# TDWG Html and ESIP or MMISW SPARQL Enpoints 
#
# T. Wellman BCB, U.S. Geological Survey 6-14-2018

from SPARQLWrapper import SPARQLWrapper, JSON
import netCDF4
import logging
import re
import os
import sys
import glob
import json
import getopt
import rdflib
import collections
from urllib import request
from datetime import datetime
try: 
    from BeautifulSoup import BeautifulSoup, Comment
except ImportError:
    from bs4 import BeautifulSoup, Comment


# In[2]:


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


# In[3]:


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
    # print("\n\t{}{}\n\t{}{}".format('TDWG vocab terms: ', orig_dic_length, 'ESIP vocab terms: ', add_terms))
    
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


# In[4]:


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
            if o in '--erddapdir':
                commands['erddap_data_dir'] = a
            elif o in '--proc_overwrite':
                commands['proc_overwrite'] = a
            elif o in '--report_dir':
                commands['report_dir'] = a
            elif o in '--help':
                get_ipython().system('usage()')
                exit(0)
                
    return commands


# In[5]:


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
                    
                logger.debug("{}{}".format('match type summary, fname: ', ncfile.split('/')[-1]))
                logger.debug("\t{}{}".format('variables - w/ existing info in file: ', srch_dic['match_count']['verbatim']))
                logger.debug("\t{}{}".format('found new var match - case insensitive: ', srch_dic['match_count']['ignore_case']))
                logger.debug("\t{}{}".format('found new var match - unique fuzzy: ', srch_dic['match_count']['fuzzy']))
                logger.debug("\t{}{}".format('found non-unique fuzzy var matches: ', srch_dic['match_count']['nonunique']))
                logger.debug("\t{}{}".format('no matches found: ',srch_dic['match_count']['no_match']))
                
        
    
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


# In[6]:


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
    # folder with netCDF files (converted from source csv)
    #
    erddap_data_dir = './erddap_data/nc_test' # nc_store'

    #
    # folder to store error and other report files 
    #
    report_dir = erddap_data_dir + '/processing_reports'

    #
    # flag whether to overwrite converted files.  If set to False, reprocessing input file will be skipped.
    #
    proc_overwrite = False
    
    
    # assemble commands dictionary - (arguments and options) 
    # ------------------------------------------------------

    commands = dict([
    ('erddap_data_dir', erddap_data_dir),
    ('report_dir', report_dir), 
    ('proc_overwrite', proc_overwrite),
    ])
    
    return commands


# In[7]:


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


# In[9]:




logform(True, logging.DEBUG) 

if __name__ == "__main__":
    
    opt_fields = ['erddapdir=', 'proc_overwrite', 'report_dir=']

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'fc', opt_fields)
    except getopt.GetoptError as err:
        logger.error(str(err))
        get_ipython().system('usage(); exit(2) ')

    # retrieve processing commands
    # -----------------------------
    logger.message("{}".format('*** Running NetCDF variable metadata checker algorithm ***'))
    commands = default_inputs() 
    commands = arg_overwrite(opts, args, commands)
    eval_vocab(commands)
    logger.message("{}".format('complete'))
    

