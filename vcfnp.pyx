# cython: profile = True

"""
Utility functions to load data from a VCF file into a numpy array.

"""

import numpy as np
cimport numpy as np
from vcflib import TYPE_FLOAT, TYPE_INTEGER, TYPE_STRING, TYPE_BOOL, TYPE_UNKNOWN
from vcflib cimport PyVariantCallFile, VariantCallFile, Variant
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from cython.operator cimport dereference as deref
import sys
import time


cdef extern from "split.h":
    # split a string on a single delimiter character (delim)
    vector[string]& split(const string &s, char delim, vector[string] &elems)
    vector[string]  split(const string &s, char delim)
    # split a string on any character found in the string of delimiters (delims)
    vector[string]& split(const string &s, const string& delims, vector[string] &elems)
    vector[string]  split(const string &s, const string& delims)


VARIANT_FIELDS = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 
                  'num_alleles', 'is_snp')


DEFAULT_DTYPE = {'CHROM': 'a12',
                 'POS': 'i4',
                 'ID': 'a12',
                 'REF': 'a12',
                 'ALT': 'a12',
                 'QUAL': 'f4',
                 'num_alleles': 'u1',
                 'is_snp': 'b1'}


DEFAULT_TYPE_MAP = {TYPE_FLOAT: 'f4',
                    TYPE_INTEGER: 'i4',
                    TYPE_STRING: 'a12',
                    TYPE_BOOL: 'b1',
                    TYPE_UNKNOWN: 'a12' # leave as string
                    }


DEFAULT_FILL_MAP = {TYPE_FLOAT: 0.,
                    TYPE_INTEGER: 0,
                    TYPE_STRING: '',
                    TYPE_BOOL: False,
                    TYPE_UNKNOWN: '' 
                    }


DEFAULT_ARITY = {'CHROM': 1,
                 'POS': 1,
                 'ID': 1,
                 'REF': 1,
                 'ALT': 2,
                 'QUAL': 1,
                 'num_alleles': 1,
                 'is_snp': 1}


DEFAULT_FILL = {'CHROM': '',
                'POS': 0,
                'ID': '',
                'REF': '',
                'ALT': '',
                'QUAL': 0,
                'num_alleles': 0,
                'is_snp': False}


cdef char SEMICOLON = ';'
cdef string FIELD_CHROM = 'CHROM'
cdef string FIELD_POS = 'POS'
cdef string FIELD_ID = 'ID'
cdef string FIELD_REF = 'REF'
cdef string FIELD_ALT = 'ALT'
cdef string FIELD_QUAL = 'QUAL'
cdef string FIELD_FILTER = 'FILTER'
cdef string FIELD_INFO = 'INFO'
cdef string FIELD_NUM_ALLELES = 'num_alleles'
cdef string FIELD_IS_SNP = 'is_snp'


def variants(filename,                  # name of VCF file
             region=None,               # region to extract
             fields=None,               # fields to extract
             dtypes=None,               # override default dtypes
             arities=None,              # override how many values to expect
             fills=None,                # override default fill values
             count=None,                # attempt to extract exactly this many records
             progress=0,                # if >0 log progress
             logstream=sys.stderr       # stream for logging progress
             ):
    
    # determine fields to extract
    if fields is None:
        fields = VARIANT_FIELDS
    else:
        for f in fields:
            assert f in VARIANT_FIELDS, 'unknown field: %s' % f
    
    # determine a numpy dtype for each field
    if dtypes is None:
        dtypes = dict()
    for f in fields:
        if f == 'FILTER':
            filterIds = PyVariantCallFile(filename).filterIds
            t = [('PASS', 'b1')]
            t += [(flt, 'b1') for flt in sorted(filterIds)]
            dtypes[f] = t
        elif f not in dtypes:
            dtypes[f] = DEFAULT_DTYPE[f]
            
    # determine expected number of values for each field
    if arities is None:
        arities = dict()
    for f in fields:
        if f == 'FILTER':
            arities[f] = 1 # one structured value
        elif f not in arities:
            arities[f] = DEFAULT_ARITY[f]
    
    # determine fill values to use where number of values is less than expectation
    if fills is None:
        fills = dict()
    for f in fields:
        if f == 'FILTER':
            fills[f] = False
        elif f not in fills:
            fills[f] = DEFAULT_FILL[f]

    # construct a numpy dtype for structured array
    dtype = list()
    for f in fields:
        t = dtypes[f]
        n = arities[f]
        if n == 1:
            dtype.append((f, t))
        else:
            dtype.append((f, t, (n,)))
            
    # set up iterator
    it = itervariants(filename, region, fields, arities, fills, progress, logstream)
    
    # build an array from the iterator
    if count is not None:
        a = np.fromiter(it, dtype=dtype, count=count)
    else:
        a = np.fromiter(it, dtype=dtype)
    return a


def itervariants(filename, 
                 region,
                 fields, 
                 arities,
                 fills,
                 progress, 
                 logstream):
    cdef VariantCallFile *variantFile
    cdef Variant *var
    
    variantFile = new VariantCallFile()
    variantFile.open(filename)
    variantFile.parseSamples = False
    if region is not None:
        variantFile.setRegion(region)
    var = new Variant(deref(variantFile))
    i = 0
    filterIds = <list>variantFile.filterIds()
    filterIds = ['PASS'] + filterIds

    if progress > 0:
        before_all = time.time()
        before = before_all
        
    while variantFile.getNextVariant(deref(var)):
        
        out = [_mkvval(var, f, arities[f], fills[f], filterIds) for f in fields]
        yield tuple(out)
        i += 1
        
        if progress > 0 and i > 0 and i % progress == 0:
            after = time.time()
            print >>logstream, '%s rows in %.2fs; batch in %.2fs (%d rows/s)' % (i, after-before_all, after-before, progress/(after-before))
            before = after

    if progress > 0:
        after_all = time.time()
        print >>logstream, '%s rows in %.2fs (%d rows/s)' % (i, after_all-before_all, i/(after_all-before_all))
        
    del variantFile
    del var

   
cdef object _mkvval(Variant *var, string field, int arity, object fill, filterIds):
    if field == FIELD_CHROM:
        out = var.sequenceName
    elif field == FIELD_POS:
        out = var.position
    elif field == FIELD_ID:
        out = var.id
    elif field == FIELD_REF:
        out = var.ref
    elif field == FIELD_ALT:
        out = _mkaltval(var, arity, fill)
    elif field == FIELD_QUAL:
        out = var.quality
    elif field == FIELD_FILTER:
        out = _mkfilterval(var, filterIds)
    elif field == FIELD_NUM_ALLELES:
        out = var.alt.size() + 1
    elif field == FIELD_IS_SNP:
        out = _is_snp(var)
    else:
        out = 0 # TODO review this
    return out
 
 
cdef object _mkaltval(Variant *var, int arity, object fill):
    if arity == 1:
        if var.alt.size() == 0:
            out = fill
        else:
            out = var.alt[0]
    elif var.alt.size() == arity:
        out = var.alt
        out = tuple(out)
    elif var.alt.size() > arity:
        out = var.alt
        out = tuple(out[:arity])
    else:
        out = var.alt
        out += [fill] * (arity-var.alt.size())
        out = tuple(out)
    return out
 
 
cdef object _mkfilterval(Variant *var, filterIds):
    filters = <list>split(var.filter, SEMICOLON)
    out = [(id in filters) for id in filterIds]
    out = tuple(out)
    return out


cdef object _is_snp(Variant *var):
    cdef int i
    cdef bytes alt
    if var.ref.size() > 1:
        return False
    for i in range(var.alt.size()):
        alt = var.alt[i]
        if alt not in {'A', 'C', 'G', 'T'}:
            return False
    return True
    

def info(filename,                  # name of VCF file
         region=None,               # region to extract
         fields=None,               # INFO fields to extract
         dtypes=None,               # override default dtypes
         arities=None,              # override how many values to expect
         fills=None,                # override default fill values
         count=None,                # attempt to extract exactly this many records
         progress=0,                # if >0 log progress
         logstream=sys.stderr       # stream for logging progress
         ):
    
    vcf = PyVariantCallFile(filename)
    infoIds = vcf.infoIds
    infoTypes = vcf.infoTypes
    infoCounts = vcf.infoCounts

    # determine INFO fields to extract
    if fields is None:
        fields = infoIds # extract all INFO fields
    else:
        for f in fields:
            assert f in infoIds, 'unknown field: %s' % f
    
    # determine a numpy dtype for each field
    if dtypes is None:
        dtypes = dict()
    for f in fields:
        if f not in dtypes:
            vcf_type = infoTypes[f]
            dtypes[f] = DEFAULT_TYPE_MAP[vcf_type]
            
    # determine expected number of values for each field
    if arities is None:
        arities = dict()
    for f in fields:
        if f not in arities:
            arities[f] = infoCounts[f]
    
    # determine fill values to use where number of values is less than expectation
    if fills is None:
        fills = dict()
    for f in fields:
        if f not in fills:
            vcf_type = infoTypes[f]
            fills[f] = DEFAULT_FILL_MAP[vcf_type]

    # construct a numpy dtype for structured array
    dtype = list()
    for f in fields:
        t = dtypes[f]
        n = arities[f]
        if n == 1:
            dtype.append((f, t))
        else:
            dtype.append((f, t, (n,)))
            
    # set up iterator
    it = iterinfo(filename, region, fields, arities, fills, progress, logstream)
    
    # build an array from the iterator
    if count is not None:
        a = np.fromiter(it, dtype=dtype, count=count)
    else:
        a = np.fromiter(it, dtype=dtype)
    return a


def iterinfo(filename, 
             region,
             fields, 
             arities,
             fills,
             progress, 
             logstream):
    cdef VariantCallFile *variantFile
    cdef Variant *var

    # TODO
    


