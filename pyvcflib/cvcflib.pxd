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
        map[string, vector[string]] info
        map[string, bool] infoFlags
        vector[string] format
        map[string, map[string, vector[string]]] samples
        vector[string] sampleNames
        vector[string] outputSampleNames
        float getAAF()
        float getNucleotideDiversity()

    cdef cppclass VariantCallFile:
        VariantCallFile()
        bint openFile(string)
        bint getNextVariant(Variant)