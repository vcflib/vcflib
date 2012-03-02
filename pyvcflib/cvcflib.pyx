from cpython cimport bool
from libcpp.vector cimport vector
from libcpp.map cimport map
from cython.operator cimport dereference as deref

cdef extern from *:
    ctypedef char* const_char_ptr "const char*"

cdef extern from "<string>" namespace "std":
    cdef cppclass string:
        string()
        string(char *)
        #char *c_str()
        const_char_ptr c_str()

        bint operator==(string&)
        bint operator==(char*)

cdef extern from "Variant.h" namespace "vcf":
    cdef cppclass Variant:
        Variant()
        unsigned long position
        string id
        string ref
        vector[string] alt
        vector[string] alleles
        string filter
        double quality
        float getAAF()
        float getNucleotideDiversity()

    cdef cppclass VariantCallFile:
        VariantCallFile()
        bint openFile(string)
        bint getNextVariant(Variant)


cdef class PyVariant:
    cdef Variant *thisptr
    def __init__(self, v):
        self.thisptr = v
    def __dealloc__(self):
        del self.thisptr

    property chrom:
        """ the chromosome of the variant"""
        def __get__(self):
            return self.thisptr.ref.c_str()

cdef PyVariant create_variant(Variant *v):
    cdef PyVariant pyvar = PyVariant.__new__(PyVariant)
    pyvar.thisptr = v
    return pyvar

cdef class VariantFile:
    cdef VariantCallFile *vcffile_ptr
    
    def __cinit__(self, vcf_file):
        self.vcffile_ptr = new VariantCallFile()
        self.vcffile_ptr.openFile(string(vcf_file))
    
    def __dealloc__(self):
        del self.vcffile_ptr
        
    def __iter__(self):
        return self

    def __next__(self):
        cdef Variant *variant = new Variant()
        success = self.vcffile_ptr.getNextVariant(deref(variant))
        if success:
            return create_variant(variant)
        else:
            raise StopIteration
