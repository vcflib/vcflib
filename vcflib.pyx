"""
Cython wrapper for classes defined in Variant.cpp.

Try to keep it simple and stay close to the C++ API.

"""


from cython.operator cimport dereference as deref
from collections import namedtuple


VariantTuple = namedtuple('Variant', ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'samples'])


# expose constants to Python
TYPE_FLOAT = FIELD_FLOAT
TYPE_INTEGER = FIELD_INTEGER
TYPE_BOOL = FIELD_BOOL
TYPE_STRING = FIELD_STRING
TYPE_UNKNOWN = FIELD_UNKNOWN


cdef class PyVariantCallFile:

    cdef VariantCallFile *thisptr

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
        var = new Variant(deref(self.thisptr))
        while self.thisptr.getNextVariant(deref(var)):
            yield VariantTuple(var.sequenceName, 
                               var.position, 
                               var.id, 
                               var.ref, 
                               var.alt, 
                               var.quality, 
                               var.filter,
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
        
