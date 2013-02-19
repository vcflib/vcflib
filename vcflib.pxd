from libcpp cimport bool
from libcpp.string cimport string


cdef extern from "Variant.h" namespace "vcf":

    cdef cppclass VariantCallFile:
        bool open(string& filename)
        bool is_open()
        bool getNextVariant(Variant& var)
        string header

    cdef cppclass Variant:
        Variant()
        Variant(VariantCallFile& v)
        void setVariantCallFile(VariantCallFile& v)


