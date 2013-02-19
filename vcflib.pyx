#from libcpp cimport bool
#from libcpp.string cimport string
#
#
#cdef extern from "Variant.h" namespace "vcf":
#
#    cdef cppclass VariantCallFile:
#        bool open(string& filename)
#        bool is_open()
#        bool getNextVariant(Variant& var)
#        string header
#
#    cdef cppclass Variant:
#        Variant()
#        Variant(VariantCallFile& v)
#        void setVariantCallFile(VariantCallFile& v)


def count_variants(filename):
    cdef VariantCallFile variantFile
    cdef string fn = filename
    cdef Variant var
    variantFile.open(fn);
    if not variantFile.is_open():
        raise Exception('variant call file is not open')
    var.setVariantCallFile(variantFile)
    n = 0
    while variantFile.getNextVariant(var):
        n += 1
    return n


