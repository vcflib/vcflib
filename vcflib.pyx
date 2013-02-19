
from cython.operator cimport dereference as deref


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
            yield (var.sequenceName, var.position, var.id, var.ref, var.alt, var.quality, var.filter)
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
            return self.thisptr.infoTypes

    property formatCounts:
        def __get__(self):
            return self.thisptr.formatTypes
        
    property parseSamples:
        def __get__(self):
            return self.thisptr.parseSamples
        def __set__(self, v):
            self.thisptr.parseSamples = v


