# cython: profile = True

"""
Utility functions to extract data from a VCF file and load into a numpy array.

"""


import numpy as np
cimport numpy as np
#from vcflib import TYPE_FLOAT, TYPE_INTEGER, TYPE_STRING, TYPE_BOOL, TYPE_UNKNOWN
from vcflib cimport PyVariantCallFile, VariantCallFile, Variant, VariantFieldType, FIELD_FLOAT, FIELD_INTEGER, FIELD_STRING, FIELD_BOOL, FIELD_UNKNOWN, ALLELE_NUMBER, GENOTYPE_NUMBER
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libc.stdlib cimport atoi, atof
from cython.operator cimport dereference as deref
import sys
import time
#from cython.view cimport array as cvarray

cdef size_t npos = -1

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

DEFAULT_VARIANT_DTYPE = {
                         'CHROM': 'a12',
                         'POS': 'i4',
                         'ID': 'a12',
                         'REF': 'a12',
                         'ALT': 'a12',
                         'QUAL': 'f4',
                         'num_alleles': 'u1',
                         'is_snp': 'b1',
                         }

DEFAULT_VARIANT_ARITY = {
                         'CHROM': 1,
                         'POS': 1,
                         'ID': 1,
                         'REF': 1,
                         'ALT': 1, # default assume biallelic (1 alt allele)
                         'QUAL': 1,
                         'num_alleles': 1,
                         'is_snp': 1
                         }

DEFAULT_VARIANT_FILL = {'CHROM': '',
                        'POS': 0,
                        'ID': '',
                        'REF': '',
                        'ALT': '',
                        'QUAL': 0,
                        'num_alleles': 0,
                        'is_snp': False}

DEFAULT_TYPE_MAP = {FIELD_FLOAT: 'f4',
                    FIELD_INTEGER: 'i4',
                    FIELD_STRING: 'a12',
                    FIELD_BOOL: 'b1',
                    FIELD_UNKNOWN: 'a12' # leave as string
                    }

DEFAULT_FILL_MAP = {FIELD_FLOAT: 0.,
                    FIELD_INTEGER: 0,
                    FIELD_STRING: '.',
                    FIELD_BOOL: False,
                    FIELD_UNKNOWN: '' 
                    }

# set some lower precision defaults for known INFO fields
DEFAULT_INFO_DTYPE = {
                     'AC': 'u2',
                     'AN': 'u2',
                     'HRun': 'u2',
                     'MLEAC': 'u2',
                     'MQ': 'f2',
                     'QD': 'f2',
                     'RPA': 'u2',
                     }

SAMPLE_FIELDS = ('is_called', 'is_phased', 'genotype')

DEFAULT_SAMPLE_DTYPE = {
                        'is_called': 'b1',
                        'is_phased': 'b1',
                        'genotype': 'i1',
                        # set some lower precision defaults for known FORMAT fields
                        'AD': 'u2',
                        'DP': 'u2',
                        'GQ': 'u1',
                        'MLPSAC': 'u1',
                        'MLPSAF': 'f2',
                        'MQ0': 'u2',
                        'PL': 'u2',
                        }

DEFAULT_SAMPLE_FILL = {
                       'is_called': False,
                       'is_phased': False,
                       'genotype': -1,
                       }

DEFAULT_SAMPLE_ARITY = {
                       'is_called': 1,
                       'is_phased': 1,
                       # N.B., set genotype arity to ploidy
                       }


cdef char SEMICOLON = ';'
cdef char DOT = '.'
cdef string GT_DELIMS = '/|'
cdef string FIELD_NAME_CHROM = 'CHROM'
cdef string FIELD_NAME_POS = 'POS'
cdef string FIELD_NAME_ID = 'ID'
cdef string FIELD_NAME_REF = 'REF'
cdef string FIELD_NAME_ALT = 'ALT'
cdef string FIELD_NAME_QUAL = 'QUAL'
cdef string FIELD_NAME_FILTER = 'FILTER'
cdef string FIELD_NAME_INFO = 'INFO'
cdef string FIELD_NAME_NUM_ALLELES = 'num_alleles'
cdef string FIELD_NAME_IS_SNP = 'is_snp'
cdef string FIELD_NAME_IS_CALLED = 'is_called'
cdef string FIELD_NAME_IS_PHASED = 'is_phased'
cdef string FIELD_NAME_GENOTYPE = 'genotype'
cdef string FIELD_NAME_GT = 'GT'



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
    """
    Load an numpy structured array with data from the fixed fields of a VCF file 
    (excluding INFO). E.g.::
    
        >>> from vcfnp import variants
        >>> a = variants('sample.vcf')
        >>> a
        array([ ('19', 111, '.', 'A', 'C', 9.600000381469727, (False, False, False), 2, True),
               ('19', 112, '.', 'A', 'G', 10.0, (False, False, False), 2, True),
               ('20', 14370, 'rs6054257', 'G', 'A', 29.0, (True, False, False), 2, True),
               ('20', 17330, '.', 'T', 'A', 3.0, (False, True, False), 2, True),
               ('20', 1110696, 'rs6040355', 'A', 'G', 67.0, (True, False, False), 3, True),
               ('20', 1230237, '.', 'T', '.', 47.0, (True, False, False), 2, False),
               ('20', 1234567, 'microsat1', 'G', 'GA', 50.0, (True, False, False), 3, False),
               ('20', 1235237, '.', 'T', '.', 0.0, (False, False, False), 2, False),
               ('X', 10, 'rsTest', 'AC', 'A', 10.0, (True, False, False), 3, False)], 
              dtype=[('CHROM', '|S12'), ('POS', '<i4'), ('ID', '|S12'), ('REF', '|S12'), ('ALT', '|S12'), ('QUAL', '<f4'), ('FILTER', [('PASS', '|b1'), ('q10', '|b1'), ('s50', '|b1')]), ('num_alleles', '|u1'), ('is_snp', '|b1')])
        >>> a['QUAL']
        array([  9.60000038,  10.        ,  29.        ,   3.        ,
                67.        ,  47.        ,  50.        ,   0.        ,  10.        ], dtype=float32)
        >>> a['FILTER']['PASS']
        array([False, False,  True, False,  True,  True,  True, False,  True], dtype=bool)

    """
    
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
            dtypes[f] = DEFAULT_VARIANT_DTYPE[f]
            
    # determine expected number of values for each field
    if arities is None:
        arities = dict()
    for f in fields:
        if f == 'FILTER':
            arities[f] = 1 # one structured value
        elif f not in arities:
            arities[f] = DEFAULT_VARIANT_ARITY[f]
    
    # determine fill values to use where number of values is less than expectation
    if fills is None:
        fills = dict()
    for f in fields:
        if f == 'FILTER':
            fills[f] = False
        elif f not in fills:
            fills[f] = DEFAULT_VARIANT_FILL[f]

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
    out = tuple([_mkvval(var, f, arities[f], fills[f], filterIds) for f in fields])
    return out


   
cdef inline object _mkvval(Variant *var, string field, int arity, object fill, list filterIds):
    if field == FIELD_NAME_CHROM:
        out = var.sequenceName
    elif field == FIELD_NAME_POS:
        out = var.position
    elif field == FIELD_NAME_ID:
        out = var.id
    elif field == FIELD_NAME_REF:
        out = var.ref
    elif field == FIELD_NAME_ALT:
        out = _mkaltval(var, arity, fill)
    elif field == FIELD_NAME_QUAL:
        out = var.quality
    elif field == FIELD_NAME_FILTER:
        out = _mkfilterval(var, filterIds)
    elif field == FIELD_NAME_NUM_ALLELES:
        out = var.alt.size() + 1
    elif field == FIELD_NAME_IS_SNP:
        out = _is_snp(var)
    else:
        out = 0 # TODO review this
    return out
 

 
cdef inline object _mkaltval(Variant *var, int arity, object fill):
    if arity == 1:
        if var.alt.size() == 0:
            out = fill
        else:
            out = var.alt.at(0)
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
        alt = var.alt.at(i)
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
    """
    Load a numpy structured array with data from the INFO field of a VCF file. 
    E.g.::
    
        >>> from vcfnp import info
        >>> a = info('sample.vcf')
        >>> a
        array([(0, 0, 0, 0, 0.0, '.', False, False),
               (0, 0, 0, 0, 0.0, '.', False, False),
               (3, 0, 0, 14, 0.5, '.', True, True),
               (3, 0, 0, 11, 0.017000000923871994, '.', False, False),
               (2, 0, 0, 10, 0.3330000042915344, 'T', True, False),
               (3, 0, 0, 13, 0.0, 'T', False, False),
               (3, 6, 3, 9, 0.0, 'G', False, False),
               (0, 0, 0, 0, 0.0, '.', False, False),
               (0, 0, 0, 0, 0.0, '.', False, False)], 
              dtype=[('NS', '<i4'), ('AN', '<u2'), ('AC', '<u2'), ('DP', '<i4'), ('AF', '<f4'), ('AA', '|S12'), ('DB', '|b1'), ('H2', '|b1')])
    
    """
    
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
            if f in DEFAULT_INFO_DTYPE:
                # known INFO field
                dtypes[f] = DEFAULT_INFO_DTYPE[f]
            else:
                vcf_type = infoTypes[f]
                dtypes[f] = DEFAULT_TYPE_MAP[vcf_type]
            
    # determine expected number of values for each field
    if arities is None:
        arities = dict()
    for f in fields:
        if f not in arities:
            vcf_count = infoCounts[f]
            if vcf_count == ALLELE_NUMBER:
                # default to 1 (biallelic)
                arities[f] = 1
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
    else:
        out = _mkval(var.info[field], arity, fill, vcf_type)
    return out



cdef inline object _mkval(vector[string]& string_vals, int arity, object fill, VariantFieldType vcf_type):
    if vcf_type == FIELD_FLOAT:
        out = _mkval_float(string_vals, arity, fill)
    elif vcf_type == FIELD_INTEGER:
        out = _mkval_int(string_vals, arity, fill)
    else:
        # make strings by default
        out = _mkval_string(string_vals, arity, fill)
    return out
 

 
cdef inline object _mkval_string(vector[string]& string_vals, int arity, string fill):
    if arity == 1:
        if string_vals.size() > 0:
            return string_vals.at(0)
        else:
            return fill
    else:
        return _mkval_string_multi(string_vals, arity, fill)



cdef inline vector[string] _mkval_string_multi(vector[string]& string_vals, int arity, string fill):
    cdef int i
    cdef vector[string] v
    for i in range(arity):
        if i < string_vals.size():
            v.push_back(string_vals.at(i))
        else:
            v.push_back(fill)
    return v



cdef inline object _mkval_float(vector[string]& string_vals, int arity, float fill):
    if arity == 1:
        out = _mkval_float_single(string_vals, fill)
    else:
        out = _mkval_float_multi(string_vals, arity, fill)
    return out



cdef inline float _mkval_float_single(vector[string]& string_vals, float fill):
    cdef float v
    if string_vals.size() > 0:
        return atof(string_vals.at(0).c_str())
    return fill
#    cdef float v = fill
#    if string_vals.size() > 0:
#        convert(string_vals.at(0), v)
#    return v



cdef inline vector[float] _mkval_float_multi(vector[string]& string_vals, int arity, float fill):
    cdef int i
    cdef vector[float] out
    for i in range(arity):
        if i < string_vals.size():
            out.push_back(atof(string_vals.at(i).c_str()))
        else:
            out.push_back(fill)
    return out
#cdef inline vector[float] _mkval_float_multi(vector[string]& string_vals, int arity, float fill):
#    cdef int i
#    cdef float v
#    cdef vector[float] out
#    for i in range(arity):
#        v = fill
#        if i < string_vals.size():
#            convert(string_vals.at(i), v)
#        out.push_back(v)
#    return out


        
cdef inline object _mkval_int(vector[string]& string_vals, int arity, int fill):
    if arity == 1:
        out = _mkval_int_single(string_vals, fill)
    else:
        out = _mkval_int_multi(string_vals, arity, fill)
    return out



cdef inline int _mkval_int_single(vector[string]& string_vals, int fill):
    cdef int v
    if string_vals.size() > 0:
        return atoi(string_vals.at(0).c_str())
    return fill
#    cdef int v = fill
#    if string_vals.size() > 0:
#        convert(string_vals.at(0), v)
#    return v



cdef inline vector[int] _mkval_int_multi(vector[string]& string_vals, int arity, int fill):
    cdef int i
    cdef vector[int] out
    for i in range(arity):
        if i < string_vals.size():
            out.push_back(atoi(string_vals.at(i).c_str()))
        else:
            out.push_back(fill)
    return out
#cdef inline vector[int] _mkval_int_multi(vector[string]& string_vals, int arity, int fill):
#    cdef int i
#    cdef int v
#    cdef vector[int] out
#    for i in range(arity):
#        v = fill
#        if i < string_vals.size():
#            convert(string_vals.at(i), v)
#        out.push_back(v)
#    return out


      
def samples(filename,                  # name of VCF file
            region=None,               # region to extract
            samples=None,              # specify which samples to extract (default all)
            ploidy=2,                  # ploidy to assume
            fields=None,               # fields to extract
            dtypes=None,               # override default dtypes
            arities=None,              # override how many values to expect
            fills=None,                # override default fill values
            count=None,                # attempt to extract exactly this many records
            progress=0,                # if >0 log progress
            logstream=sys.stderr       # stream for logging progress
            ):
    """
    Load a numpy structured array with data from the sample columns of a VCF
    file. E.g.::
    
        >>> from vcfnp import samples
        >>> a = samples('sample.vcf')
        >>> a
        array([ ((True, True, [0, 0], '0|0', 0, 0, [10, 10]), (True, True, [0, 0], '0|0', 0, 0, [10, 10]), (True, False, [0, 1], '0/1', 0, 0, [3, 3])),
               ((True, True, [0, 0], '0|0', 0, 0, [10, 10]), (True, True, [0, 0], '0|0', 0, 0, [10, 10]), (True, False, [0, 1], '0/1', 0, 0, [3, 3])),
               ((True, True, [0, 0], '0|0', 48, 1, [51, 51]), (True, True, [1, 0], '1|0', 48, 8, [51, 51]), (True, False, [1, 1], '1/1', 43, 5, [0, 0])),
               ((True, True, [0, 0], '0|0', 49, 3, [58, 50]), (True, True, [0, 1], '0|1', 3, 5, [65, 3]), (True, False, [0, 0], '0/0', 41, 3, [0, 0])),
               ((True, True, [1, 2], '1|2', 21, 6, [23, 27]), (True, True, [2, 1], '2|1', 2, 0, [18, 2]), (True, False, [2, 2], '2/2', 35, 4, [0, 0])),
               ((True, True, [0, 0], '0|0', 54, 0, [56, 60]), (True, True, [0, 0], '0|0', 48, 4, [51, 51]), (True, False, [0, 0], '0/0', 61, 2, [0, 0])),
               ((True, False, [0, 1], '0/1', 0, 4, [0, 0]), (True, False, [0, 2], '0/2', 17, 2, [0, 0]), (True, False, [1, 1], '1/1', 40, 3, [0, 0])),
               ((True, False, [0, 0], '0/0', 0, 0, [0, 0]), (True, True, [0, 0], '0|0', 0, 0, [0, 0]), (False, False, [-1, -1], './.', 0, 0, [0, 0])),
               ((True, False, [0, -1], '0', 0, 0, [0, 0]), (True, False, [0, 1], '0/1', 0, 0, [0, 0]), (True, True, [0, 2], '0|2', 0, 0, [0, 0]))], 
              dtype=[('NA00001', [('is_called', '|b1'), ('is_phased', '|b1'), ('genotype', '|i1', (2,)), ('GT', '|S3'), ('GQ', '|u1'), ('DP', '<u2'), ('HQ', '<i4', (2,))]), ('NA00002', [('is_called', '|b1'), ('is_phased', '|b1'), ('genotype', '|i1', (2,)), ('GT', '|S3'), ('GQ', '|u1'), ('DP', '<u2'), ('HQ', '<i4', (2,))]), ('NA00003', [('is_called', '|b1'), ('is_phased', '|b1'), ('genotype', '|i1', (2,)), ('GT', '|S3'), ('GQ', '|u1'), ('DP', '<u2'), ('HQ', '<i4', (2,))])])
        >>> a['NA00001']
        array([(True, True, [0, 0], '0|0', 0, 0, [10, 10]),
               (True, True, [0, 0], '0|0', 0, 0, [10, 10]),
               (True, True, [0, 0], '0|0', 48, 1, [51, 51]),
               (True, True, [0, 0], '0|0', 49, 3, [58, 50]),
               (True, True, [1, 2], '1|2', 21, 6, [23, 27]),
               (True, True, [0, 0], '0|0', 54, 0, [56, 60]),
               (True, False, [0, 1], '0/1', 0, 4, [0, 0]),
               (True, False, [0, 0], '0/0', 0, 0, [0, 0]),
               (True, False, [0, -1], '0', 0, 0, [0, 0])], 
              dtype=[('is_called', '|b1'), ('is_phased', '|b1'), ('genotype', '|i1', (2,)), ('GT', '|S3'), ('GQ', '|u1'), ('DP', '<u2'), ('HQ', '<i4', (2,))])
    
    """
    
    vcf = PyVariantCallFile(filename)
    formatIds = vcf.formatIds
    formatTypes = vcf.formatTypes
    formatCounts = vcf.formatCounts
    all_samples = vcf.sampleNames

    if samples is None:
        samples = all_samples
    else:
        for s in samples:
            assert s in all_samples, 'unknown sample: %s' % s
            
    # determine fields to extract
    if fields is None:
        fields = list(SAMPLE_FIELDS) + formatIds
    else:
        for f in fields:
            assert f in VARIANT_FIELDS or f in formatIds, 'unknown field: %s' % f
    
    # determine a numpy dtype for each field
    if dtypes is None:
        dtypes = dict()
    for f in fields:
        if f not in dtypes:
            if f == 'GT':
                dtypes[f] = 'a%d' % ((ploidy*2)-1)
            elif f in DEFAULT_SAMPLE_DTYPE:
                # known field
                dtypes[f] = DEFAULT_SAMPLE_DTYPE[f]
            else:
                vcf_type = formatTypes[f]
                dtypes[f] = DEFAULT_TYPE_MAP[vcf_type]
            
    # determine expected number of values for each field
    if arities is None:
        arities = dict()
    for f in fields:
        if f not in arities:
            if f == 'genotype':
                arities[f] = ploidy
            elif f in DEFAULT_SAMPLE_ARITY:
                arities[f] = DEFAULT_SAMPLE_ARITY[f]
            else:
                vcf_count = formatCounts[f]
                if vcf_count == ALLELE_NUMBER:
                    # default to 1 (biallelic)
                    arities[f] = 1
                elif vcf_count == GENOTYPE_NUMBER:
                    # arity = (n + p - 1) choose p (n is number of alleles; p is ploidy)
                    # default to biallelic (n = 2)
                    arities[f] = ploidy + 1
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
            if f == 'GT':
                fills[f] = '/'.join(['.'] * ploidy)
            elif f in DEFAULT_SAMPLE_FILL:
                fills[f] = DEFAULT_SAMPLE_FILL[f]
            else:
                vcf_type = formatTypes[f]
                fills[f] = DEFAULT_FILL_MAP[vcf_type]

    # construct a numpy dtype for structured array cells
    cell_dtype = list()
    for f in fields:
        t = dtypes[f]
        n = arities[f]
        if n == 1:
            cell_dtype.append((f, t))
        else:
            cell_dtype.append((f, t, (n,)))
    # construct a numpy dtype for structured array
    dtype = [(s, cell_dtype) for s in samples]
            
    # set up iterator
    it = _itersamples(filename, region, samples, ploidy, fields, arities, fills)
    
    # build an array from the iterator
    return _fromiter(it, dtype, count, progress, logstream)



def _itersamples(filename, 
                 region,
                 vector[string] samples,
                 int ploidy,
                 vector[string] fields, 
                 map[string, int] arities,
                 dict fills):
    cdef VariantCallFile *variantFile
    cdef Variant *var

    variantFile = new VariantCallFile()
    variantFile.open(filename)
    variantFile.parseSamples = True
    if region is not None:
        variantFile.setRegion(region)
    var = new Variant(deref(variantFile))

    while variantFile.getNextVariant(deref(var)):
        yield _mkssvals(var, samples, ploidy, fields, arities, fills, variantFile.formatTypes)
#        out = [_mksvals(var, s, ploidy, fields, arities, fills, variantFile.formatTypes) for s in samples]
#        yield tuple(out)
        
    del variantFile
    del var
    
    
cdef inline object _mkssvals(Variant *var,
                             vector[string] samples,
                             int ploidy,
                             vector[string] fields, 
                             map[string, int] arities,
                             dict fills,
                             map[string, VariantFieldType]& formatTypes):
    out = [_mksvals(var, s, ploidy, fields, arities, fills, formatTypes) for s in samples]
    return tuple(out)

    
    
cdef inline object _mksvals(Variant *var, 
                            string sample,
                            int ploidy,
                            vector[string] fields, 
                            map[string, int] arities,
                            dict fills,
                            map[string, VariantFieldType]& formatTypes):
    out = [_mksval(var.samples[sample], ploidy, f, arities[f], fills[f], formatTypes) for f in fields]
    return tuple(out)
    


cdef inline object _mksval(map[string, vector[string]]& sample_data, 
                           int ploidy,
                           string field,
                           int arity, 
                           object fill, 
                           map[string, VariantFieldType]& formatTypes):
    if field == FIELD_NAME_IS_CALLED:
        return _is_called(sample_data)
    elif field == FIELD_NAME_IS_PHASED:
        return _is_phased(sample_data)
    elif field == FIELD_NAME_GENOTYPE:
        return _genotype(sample_data, ploidy)
    else:
        return _mkval(sample_data[field], arity, fill, formatTypes[field])
    


cdef inline bool _is_called(map[string, vector[string]]& sample_data):
    cdef vector[string] *gts
    gts = &sample_data[FIELD_NAME_GT]
    if gts.size() == 0:
        return False
    else:
        return (gts.at(0).find('.') == npos)
        
        
cdef inline bool _is_phased(map[string, vector[string]]& sample_data):
    cdef vector[string] *gts
    gts = &sample_data[FIELD_NAME_GT]
    if gts.size() == 0:
        return False
    else:
        return (gts.at(0).find('|') != npos)


cdef inline object _genotype(map[string, vector[string]]& sample_data, int ploidy):
    cdef vector[string] *gts
    cdef vector[int] alleles
    cdef vector[string] allele_strings
    cdef int i
    cdef int allele
    gts = &sample_data[FIELD_NAME_GT]
    if gts.size() == 0:
        if ploidy == 1:
            return -1
        else:
            return (-1,) * ploidy
    else:
        split(gts.at(0), GT_DELIMS, allele_strings)
        if ploidy == 1:
            if allele_strings.size() > 0:
                return atoi(allele_strings.at(0).c_str())
            else:
                return -1
        else:
            for i in range(ploidy):
                if i < allele_strings.size():
                    alleles.push_back(atoi(allele_strings.at(i).c_str()))
                else:
                    alleles.push_back(-1)
            return tuple(alleles)
#cdef inline object _genotype(map[string, vector[string]]& sample_data, int ploidy):
#    cdef vector[string] *gts
#    cdef vector[int] alleles
#    cdef vector[string] allele_strings
#    cdef int i
#    cdef int allele
#    gts = &sample_data[FIELD_NAME_GT]
#    if gts.size() == 0:
#        if ploidy == 1:
#            return -1
#        else:
#            return (-1,) * ploidy
#    else:
#        split(gts.at(0), GT_DELIMS, allele_strings)
#        if ploidy == 1:
#            allele = -1
#            if allele_strings.size() > 0:
#                convert(allele_strings.at(0), allele)
#            return allele
#        else:
#            for i in range(ploidy):
#                allele = -1
#                if i < allele_strings.size():
#                    convert(allele_strings.at(i), allele)
#                alleles.push_back(allele)
#            return tuple(alleles)
        
        
def view2d(a):
    """
    Utility function to view a structured 1D array where all fields have a uniform dtype 
    (e.g., an array constructed by :func:samples) as a 2D array. E.g.::
    
        >>> from vcfnp import samples
        >>> a = samples('sample.vcf')
        >>> a
        array([ ((True, True, [0, 0], '0|0', 0, 0, [10, 10]), (True, True, [0, 0], '0|0', 0, 0, [10, 10]), (True, False, [0, 1], '0/1', 0, 0, [3, 3])),
               ((True, True, [0, 0], '0|0', 0, 0, [10, 10]), (True, True, [0, 0], '0|0', 0, 0, [10, 10]), (True, False, [0, 1], '0/1', 0, 0, [3, 3])),
               ((True, True, [0, 0], '0|0', 48, 1, [51, 51]), (True, True, [1, 0], '1|0', 48, 8, [51, 51]), (True, False, [1, 1], '1/1', 43, 5, [0, 0])),
               ((True, True, [0, 0], '0|0', 49, 3, [58, 50]), (True, True, [0, 1], '0|1', 3, 5, [65, 3]), (True, False, [0, 0], '0/0', 41, 3, [0, 0])),
               ((True, True, [1, 2], '1|2', 21, 6, [23, 27]), (True, True, [2, 1], '2|1', 2, 0, [18, 2]), (True, False, [2, 2], '2/2', 35, 4, [0, 0])),
               ((True, True, [0, 0], '0|0', 54, 0, [56, 60]), (True, True, [0, 0], '0|0', 48, 4, [51, 51]), (True, False, [0, 0], '0/0', 61, 2, [0, 0])),
               ((True, False, [0, 1], '0/1', 0, 4, [0, 0]), (True, False, [0, 2], '0/2', 17, 2, [0, 0]), (True, False, [1, 1], '1/1', 40, 3, [0, 0])),
               ((True, False, [0, 0], '0/0', 0, 0, [0, 0]), (True, True, [0, 0], '0|0', 0, 0, [0, 0]), (False, False, [-1, -1], './.', 0, 0, [0, 0])),
               ((True, False, [0, -1], '0', 0, 0, [0, 0]), (True, False, [0, 1], '0/1', 0, 0, [0, 0]), (True, True, [0, 2], '0|2', 0, 0, [0, 0]))], 
              dtype=[('NA00001', [('is_called', '|b1'), ('is_phased', '|b1'), ('genotype', '|i1', (2,)), ('GT', '|S3'), ('GQ', '|u1'), ('DP', '<u2'), ('HQ', '<i4', (2,))]), ('NA00002', [('is_called', '|b1'), ('is_phased', '|b1'), ('genotype', '|i1', (2,)), ('GT', '|S3'), ('GQ', '|u1'), ('DP', '<u2'), ('HQ', '<i4', (2,))]), ('NA00003', [('is_called', '|b1'), ('is_phased', '|b1'), ('genotype', '|i1', (2,)), ('GT', '|S3'), ('GQ', '|u1'), ('DP', '<u2'), ('HQ', '<i4', (2,))])])
        >>> from vcfnp import view2d
        >>> b = view2d(a)
        >>> b
        array([[(True, True, [0, 0], '0|0', 0, 0, [10, 10]),
                (True, True, [0, 0], '0|0', 0, 0, [10, 10]),
                (True, False, [0, 1], '0/1', 0, 0, [3, 3])],
               [(True, True, [0, 0], '0|0', 0, 0, [10, 10]),
                (True, True, [0, 0], '0|0', 0, 0, [10, 10]),
                (True, False, [0, 1], '0/1', 0, 0, [3, 3])],
               [(True, True, [0, 0], '0|0', 48, 1, [51, 51]),
                (True, True, [1, 0], '1|0', 48, 8, [51, 51]),
                (True, False, [1, 1], '1/1', 43, 5, [0, 0])],
               [(True, True, [0, 0], '0|0', 49, 3, [58, 50]),
                (True, True, [0, 1], '0|1', 3, 5, [65, 3]),
                (True, False, [0, 0], '0/0', 41, 3, [0, 0])],
               [(True, True, [1, 2], '1|2', 21, 6, [23, 27]),
                (True, True, [2, 1], '2|1', 2, 0, [18, 2]),
                (True, False, [2, 2], '2/2', 35, 4, [0, 0])],
               [(True, True, [0, 0], '0|0', 54, 0, [56, 60]),
                (True, True, [0, 0], '0|0', 48, 4, [51, 51]),
                (True, False, [0, 0], '0/0', 61, 2, [0, 0])],
               [(True, False, [0, 1], '0/1', 0, 4, [0, 0]),
                (True, False, [0, 2], '0/2', 17, 2, [0, 0]),
                (True, False, [1, 1], '1/1', 40, 3, [0, 0])],
               [(True, False, [0, 0], '0/0', 0, 0, [0, 0]),
                (True, True, [0, 0], '0|0', 0, 0, [0, 0]),
                (False, False, [-1, -1], './.', 0, 0, [0, 0])],
               [(True, False, [0, -1], '0', 0, 0, [0, 0]),
                (True, False, [0, 1], '0/1', 0, 0, [0, 0]),
                (True, True, [0, 2], '0|2', 0, 0, [0, 0])]], 
              dtype=[('is_called', '|b1'), ('is_phased', '|b1'), ('genotype', '|i1', (2,)), ('GT', '|S3'), ('GQ', '|u1'), ('DP', '<u2'), ('HQ', '<i4', (2,))])
        >>> b['GT']
        array([['0|0', '0|0', '0/1'],
               ['0|0', '0|0', '0/1'],
               ['0|0', '1|0', '1/1'],
               ['0|0', '0|1', '0/0'],
               ['1|2', '2|1', '2/2'],
               ['0|0', '0|0', '0/0'],
               ['0/1', '0/2', '1/1'],
               ['0/0', '0|0', './.'],
               ['0', '0/1', '0|2']], 
              dtype='|S3')

    """
    
    rows = a.size
    cols = len(a.dtype)
    dtype = a.dtype[0]
    b = a.view(dtype).reshape(rows, cols)
    return b
    





