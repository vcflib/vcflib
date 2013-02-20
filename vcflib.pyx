"""
Cython wrapper for classes defined in Variant.cpp.

Try to keep it simple and stay close to the C++ API.

"""


from cython.operator cimport dereference as deref
from collections import namedtuple
import numpy as np
cimport numpy as np
import time
import sys


VariantTuple = namedtuple('Variant', ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'samples'])


# expose constants to Python
TYPE_FLOAT = FIELD_FLOAT
TYPE_INTEGER = FIELD_INTEGER
TYPE_BOOL = FIELD_BOOL
TYPE_STRING = FIELD_STRING
TYPE_UNKNOWN = FIELD_UNKNOWN


cdef class PyVariantCallFile:

    def __cinit__(self, filename):
        self.thisptr = new VariantCallFile()
        self.thisptr.open(filename)

    def __dealloc__(self):
        del self.thisptr        
        
    def __len__(self):
        cdef Variant var
        var.setVariantCallFile(self.thisptr)
        n = 0
        while self.thisptr.getNextVariant(var):
            n += 1
        return n

    def __iter__(self):
        cdef Variant *var
        cdef vector[string] filters
        cdef char semicolon = ';'
        var = new Variant(deref(self.thisptr))
        while self.thisptr.getNextVariant(deref(var)):
            # split the filter field here in C++ to avoid having to do it in Python later
            filters = split(var.filter, semicolon)
            yield VariantTuple(var.sequenceName, 
                               var.position, 
                               var.id, 
                               var.ref, 
                               var.alt, 
                               var.quality, 
                               filters,
                               var.info,
                               var.samples)
        del var

    property infoIds:
        def __get__(self):
            return self.thisptr.infoIds()

    property formatIds:
        def __get__(self):
            return self.thisptr.formatIds()

    property filterIds:
        def __get__(self):
            return self.thisptr.filterIds()

    property infoTypes:
        def __get__(self):
            return self.thisptr.infoTypes

    property formatTypes:
        def __get__(self):
            return self.thisptr.formatTypes

    property infoCounts:
        def __get__(self):
            return self.thisptr.infoCounts

    property formatCounts:
        def __get__(self):
            return self.thisptr.formatCounts
        
    property parseSamples:
        def __get__(self):
            return self.thisptr.parseSamples
        def __set__(self, v):
            self.thisptr.parseSamples = v

    property header:
        def __get__(self):
            return self.thisptr.header
        
    property fileformat: # [sic] no camel case
        def __get__(self):
            return self.thisptr.fileformat
        
    property fileDate:
        def __get__(self):
            return self.thisptr.fileDate
        
    property source:
        def __get__(self):
            return self.thisptr.source
        
    property reference:
        def __get__(self):
            return self.thisptr.reference
        
    property phasing:
        def __get__(self):
            return self.thisptr.phasing
        
    property sampleNames:
        def __get__(self):
            return self.thisptr.sampleNames
        
        
def itervariants(filename, int progress=0, logstream=sys.stderr):
    cdef VariantCallFile *variantFile = new VariantCallFile()
    cdef Variant *var
    cdef int i = 0
    cdef vector[string] filters
    cdef char semicolon = ';'
    variantFile.open(filename)
    variantFile.parseSamples = False
    var = new Variant(deref(variantFile))
    if progress > 0:
        before_all = time.time()
        before = before_all
    while variantFile.getNextVariant(deref(var)):
        # split the filter field here in C++ to avoid having to do it in Python later
        filters = split(var.filter, semicolon)
        yield (var.sequenceName, var.position, var.id, var.ref, tuple(var.alt), var.quality, tuple(filters))
        if progress > 0 and i > 0 and i % progress == 0:
            after = time.time()
            print >>logstream, '%s rows in %.3fs; batch in %.3fs (%d rows/s)' % (i, after-before_all, after-before, progress/(after-before))
            before = after
        i += 1
    del variantFile
    del var
    if progress > 0:
        after_all = time.time()
        print >>logstream, '%s rows in %.3fs (%d rows/s)' % (i, after_all-before_all, i/(after_all-before_all))
    
    
