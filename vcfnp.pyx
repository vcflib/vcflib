# cython: profile = True

"""
Utility functions to load data from a VCF file into a numpy array.

"""

import numpy as np
cimport numpy as np
#from vcflib import TYPE_FLOAT, TYPE_INTEGER, TYPE_STRING, TYPE_BOOL, TYPE_UNKNOWN
from vcflib cimport PyVariantCallFile, VariantCallFile, Variant, VariantFieldType, FIELD_FLOAT, FIELD_INTEGER, FIELD_STRING, FIELD_BOOL, FIELD_UNKNOWN, ALLELE_NUMBER, GENOTYPE_NUMBER
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
    

cdef extern from "convert.h":
    bool convert(const string& s, int& r)
    bool convert(const string& s, float& r)


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


DEFAULT_TYPE_MAP = {FIELD_FLOAT: 'f4',
                    FIELD_INTEGER: 'i4',
                    FIELD_STRING: 'a12',
                    FIELD_BOOL: 'b1',
                    FIELD_UNKNOWN: 'a12' # leave as string
                    }


DEFAULT_FILL_MAP = {FIELD_FLOAT: 0.,
                    FIELD_INTEGER: 0,
                    FIELD_STRING: '',
                    FIELD_BOOL: False,
                    FIELD_UNKNOWN: '' 
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
cdef string ATTR_CHROM = 'CHROM'
cdef string ATTR_POS = 'POS'
cdef string ATTR_ID = 'ID'
cdef string ATTR_REF = 'REF'
cdef string ATTR_ALT = 'ALT'
cdef string ATTR_QUAL = 'QUAL'
cdef string ATTR_FILTER = 'FILTER'
cdef string ATTR_INFO = 'INFO'
cdef string ATTR_NUM_ALLELES = 'num_alleles'
cdef string ATTR_IS_SNP = 'is_snp'


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
    it = _itervariants(filename, region, fields, arities, fills)
    
    # build an array from the iterator
    return _fromiter(it, dtype, count, progress, logstream)


def _fromiter(it, dtype, count, int progress=0, logstream=sys.stderr):
    if progress > 0:
        it = _iter_withprogress(it, progress, logstream)
    if count is not None:
        a = np.fromiter(it, dtype=dtype, count=count)
    else:
        a = np.fromiter(it, dtype=dtype)
    return a


def _iter_withprogress(iterable, int progress, logstream):
    cdef int i
    before_all = time.time()
    before = before_all
    for i, o in enumerate(iterable):
        yield o
        if i > 0 and i % progress == 0:
            after = time.time()
            print >>logstream, '%s rows in %.2fs; batch in %.2fs (%d rows/s)' % (i, after-before_all, after-before, progress/(after-before))
            before = after
    after_all = time.time()
    print >>logstream, '%s rows in %.2fs (%d rows/s)' % (i, after_all-before_all, i/(after_all-before_all))


def _itervariants(filename, 
                 region,
                 vector[string] fields, 
                 map[string, int] arities,
                 dict fills):
    cdef VariantCallFile *variantFile
    cdef Variant *var
    cdef vector[string] filterIds
    
    variantFile = new VariantCallFile()
    variantFile.open(filename)
    variantFile.parseSamples = False
    if region is not None:
        variantFile.setRegion(region)
    var = new Variant(deref(variantFile))
    filterIds = <list>variantFile.filterIds()
    filterIds = ['PASS'] + filterIds

    while variantFile.getNextVariant(deref(var)):
        yield _mkvvals(var, fields, arities, fills, filterIds)
        
    del variantFile
    del var


cdef inline object _mkvvals(Variant *var, 
                            vector[string] fields, 
                            map[string, int] arities, 
                            dict fills, 
                            list filterIds):
    out = [_mkvval(var, f, arities[f], fills[f], filterIds) for f in fields]
    return tuple(out)

   
cdef inline object _mkvval(Variant *var, string field, int arity, object fill, list filterIds):
    if field == ATTR_CHROM:
        out = var.sequenceName
    elif field == ATTR_POS:
        out = var.position
    elif field == ATTR_ID:
        out = var.id
    elif field == ATTR_REF:
        out = var.ref
    elif field == ATTR_ALT:
        out = _mkaltval(var, arity, fill)
    elif field == ATTR_QUAL:
        out = var.quality
    elif field == ATTR_FILTER:
        out = _mkfilterval(var, filterIds)
    elif field == ATTR_NUM_ALLELES:
        out = var.alt.size() + 1
    elif field == ATTR_IS_SNP:
        out = _is_snp(var)
    else:
        out = 0 # TODO review this
    return out
 
 
cdef inline object _mkaltval(Variant *var, int arity, object fill):
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
 
 
cdef inline object _mkfilterval(Variant *var, list filterIds):
    filters = <list>split(var.filter, SEMICOLON)
    out = [(id in filters) for id in filterIds]
    out = tuple(out)
    return out


cdef inline object _is_snp(Variant *var):
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
            vcf_count = infoCounts[f]
            if vcf_count == ALLELE_NUMBER:
                # can't deal with variable arity, default to 2 (biallelic)
                arities[f] = 2
            elif vcf_count <= 0:
                # catch any other cases of non-specific arity
                arities[f] = 1
            else:
                arities[f] = vcf_count
    
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
    it = _iterinfo(filename, region, fields, arities, fills)
    
    # build an array from the iterator
    return _fromiter(it, dtype, count, progress, logstream)


def _iterinfo(filename, 
             region,
             vector[string] fields, 
             map[string, int] arities,
             dict fills):
    cdef VariantCallFile *variantFile
    cdef Variant *var

    variantFile = new VariantCallFile()
    variantFile.open(filename)
    variantFile.parseSamples = False
    if region is not None:
        variantFile.setRegion(region)
    var = new Variant(deref(variantFile))

    while variantFile.getNextVariant(deref(var)):
        yield _mkivals(var, fields, arities, fills, variantFile.infoTypes)
        
    del variantFile
    del var
    
    
cdef inline object _mkivals(Variant *var, 
                            vector[string] fields, 
                            map[string, int] arities,
                            dict fills,
                            map[string, VariantFieldType] infoTypes):
    out = [_mkival(var, f, arities[f], fills[f], infoTypes[f]) for f in fields]
    return tuple(out)
    
    
cdef inline object _mkival(Variant *var, string field, int arity, object fill, VariantFieldType vcf_type):
    if vcf_type == FIELD_BOOL:
        # ignore arity, this is a flag
        out = (var.infoFlags.count(field) > 0)
    elif vcf_type == FIELD_STRING:
        out = _mkival_string(var, field, arity, fill)
    elif vcf_type == FIELD_FLOAT:
        out = _mkival_float(var, field, arity, fill)
    elif vcf_type == FIELD_INTEGER:
        out = _mkival_int(var, field, arity, fill)
    else:
        # fall back to strings
        out = _mkival_string(var, field, arity, fill)
    return out
 
 
cdef inline object _mkival_string(Variant *var, string field, int arity, string fill):
    if arity == 1:
        out = _mkival_string_single(&var.info[field], fill)
    else:
        out = _mkival_string_multi(&var.info[field], arity, fill)
    return out


cdef inline string _mkival_string_single(vector[string] *string_vals, string fill):
    if string_vals.size() > 0:
        return string_vals.at(0)
    else:
        return fill


cdef inline vector[string] _mkival_string_multi(vector[string] *string_vals, int arity, string fill):
    cdef int i
    cdef vector[string] v
    for i in range(arity):
        if i < string_vals.size():
            v.push_back(string_vals.at(i))
        else:
            v.push_back(fill)
    return v


cdef inline object _mkival_float(Variant *var, string field, int arity, float fill):
    if arity == 1:
        out = _mkival_float_single(&var.info[field], fill)
    else:
        out = _mkival_float_multi(&var.info[field], arity, fill)
    return out


cdef inline float _mkival_float_single(vector[string] *string_vals, float fill):
    cdef float v = fill
    if string_vals.size() > 0:
        convert(string_vals.at(0), v)
    return v


cdef inline object _mkival_float_multi(vector[string] *string_vals, int arity, float fill) except +:
    cdef int i
    cdef float v
    cdef vector[float] out
    for i in range(arity):
        v = fill
        if i < string_vals.size():
            convert(string_vals.at(i), v)
        out.push_back(v)
    return out

        
cdef inline object _mkival_int(Variant *var, string field, int arity, int fill):
    if arity == 1:
        out = _mkival_int_single(&var.info[field], fill)
    else:
        out = _mkival_int_multi(&var.info[field], arity, fill)
    return out


cdef inline int _mkival_int_single(vector[string] *string_vals, int fill):
    cdef int v = fill
    if string_vals.size() > 0:
        convert(string_vals.at(0), v)
    return v


cdef inline object _mkival_int_multi(vector[string] *string_vals, int arity, int fill) except +:
    cdef int i
    cdef int v
    cdef vector[int] out
    for i in range(arity):
        v = fill
        if i < string_vals.size():
            convert(string_vals.at(i), v)
        out.push_back(v)
    return out

        


