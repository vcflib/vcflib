from cpython cimport bool
from libcpp.vector cimport vector
from libcpp.map cimport map
from cython.operator cimport dereference as deref, preincrement as inc

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


###############################################
# Cython definitions of the C++ vcflib classes
# 1. VariantCallFile: exposed to Python as
#    VariantFile (below)
# 2. Variant: exposed to Python as
#    VariantRecord (below)
###############################################
cdef extern from "Variant.h" namespace "vcf":
    cdef cppclass VariantCallFile:
        VariantCallFile()
        bint open(string &)
        bint getNextVariant(Variant)
        string header # contains the entire header

    cdef cppclass Variant:
        Variant(VariantCallFile &)
        string echo() # "print" Variant to string. for __repr__

        # core VCF fields
        string sequenceName
        unsigned long position
        string id
        string ref
        vector[string] alt
        double quality
        string filter
        map[string, vector[string] ] info
        vector[string] format
        
        # derived fields
        vector[string] alleles
        map[string, map[string, vector[string]]] samples
        vector[string] sampleNames
        
        unsigned int num_hom_ref     # how many hom_ref gts are there?
        unsigned int num_het         # how many het gts are there?
        unsigned int num_hom_alt     # how many hom_alt gts are there?
        unsigned int num_unknown     # how many unknown gts are there?
        unsigned int num_valid       # how many valid (i.e., !unknown) gts are there?
        unsigned int num_alt_alleles # how many alternate alleles are there?
        vector[string] gts           # a vector of the nucl. genotypes for each
                                     # sample. (e.g. [AA AG GG]). in sample order
        vector[int] gt_types         # a vector of the genotype classes for each
                                     # sample. (e.g. [0 1 2). in sample order
        # TO DO.  Figure out why Cython doesn't handle bool here.
        vector[int] gt_phases        # a vector of booleans described the phase 
                                     # of each gt for each sample (in order)
        # custom methods
        float getAAF()
        float getNucleotideDiversity()



###############################################
# Utility function to convert from STL data
# structures (e.g., vector and map) to
# Python data structures (e.g., list and dict)
#
# TODO: how to template vec2list?
###############################################
cdef list string_vec2list(vector[string] sv):
    """
    convert an STL vector<string> to a list
    """
    cdef size_t size = sv.size(), i
    return [sv.at(i).c_str() for i in range(size)]


cdef list int_vec2list(vector[int] sv):
    """
    convert an STL vector<int> to a list
    """
    cdef size_t size = sv.size(), i
    return [sv.at(i) for i in range(size)]


cdef list bool_vec2list(vector[int] sv):
    """
    convert an STL vector<bool> to a list
    """
    cdef size_t size = sv.size(), i
    return [sv.at(i) > 0 for i in range(size)]


cdef dict string_map2dict(map[string, vector[string] ] sm):
    """
    convert an STL map<string, vector<string> > to a dict
    """
    dict = {}
    cdef map[string, vector[string] ].iterator it = sm.begin()
    while it != sm.end():
        dict[deref(it).first.c_str()] = string_vec2list(deref(it).second)
        inc(it)
    return dict


cdef VariantRecord create_variant(Variant *v):
    """
    A utility function to store a C++ Variant object
    in a Python object that can be accessed by the user.
    """
    cdef VariantRecord pyvar = VariantRecord.__new__(VariantRecord)
    pyvar._thisptr = v
    return pyvar


################################################
# The Python interface, which consists of the 
# following classes:
# 1. VariantFile - open and iterate through a
#       VCF file.  each iteration yields a
#       instance of VariantRecord (#2)
# 2. VariantRecord - a class containing the core
#       VCF variant attributes (e.g., REF, ALT)
#       as well as many other useful, derived 
#       attributes (e.g., aaf, pi, etc.)
################################################
cdef class VariantFile:
    """
    A Python class for iterating through each record in
    a VCF file.
    
    Example usage:
       vcf_file = pyvcflib.VariantFile("sample.vcf")
       for var in vcf_file:
           print var.chrom, var.pos, var.ref, var.alt
    """
    cdef VariantCallFile *vcffile_ptr

    def __cinit__(self, vcf_file):
        self.vcffile_ptr = new VariantCallFile()
        self.vcffile_ptr.open(string(vcf_file))

    def __dealloc__(self):
        del self.vcffile_ptr

    def __iter__(self):
        return self

    def __next__(self):
        cdef Variant *variant = new Variant(deref(self.vcffile_ptr))
        success = self.vcffile_ptr.getNextVariant(deref(variant))
        if success:
            return create_variant(variant)
        else:
            raise StopIteration

    property header:
        """ the entire VCF file's header"""
        def __get__(self):
            return self.vcffile_ptr.header.c_str()


cdef class VariantRecord:
    """
    A Python wrapper class for vcflib's Variant class.
    """
    cdef Variant *_thisptr

    def __init__(self):
        pass

    def __dealloc__(self):
        del self._thisptr

    def __repr__(self):
        o = self._thisptr.echo()
        return o.c_str()
        
    def __hash__(self):
        #TODO
        pass

    # core attributes from VCF spec
    property chrom:
        """ the chromosome of the variant"""
        def __get__(self):
            return self._thisptr.sequenceName.c_str()

    property pos:
        """ the 1-based position of the variant"""
        def __get__(self):
            return self._thisptr.position

    property id:
        """ the database id for the variant"""
        def __get__(self):
            return self._thisptr.id.c_str()

    property ref:
        """ the reference allele for the variant"""
        def __get__(self):
            return self._thisptr.ref.c_str()

    property alt:
        """ the alternate allele for the variant"""
        def __get__(self):
            return string_vec2list(self._thisptr.alt)

    property qual:
        """ the PHRED-like variant quality estimate"""
        def __get__(self):
            return self._thisptr.quality

    property filter:
        """ the filters that are applicable to this variant"""
        def __get__(self):
            return self._thisptr.filter.c_str()

    property info:
        """ the set of annotations variant"""
        def __get__(self):
            return string_map2dict(self._thisptr.info)

    property format:
        """ the format of the genotypes for this variant"""
        def __get__(self):
            return string_vec2list(self._thisptr.format)

    # dervived attributes
    property sample_names:
        """ the list of sample names"""
        def __get__(self):
            return string_vec2list(self._thisptr.sampleNames)

    property samples:
        """a dictionary of genotype attributes for each sample"""
        def __get__(self):
            dict = {}
            samples = self._thisptr.samples
            cdef map[string, map[string, vector[string]]].iterator it = \
                samples.begin()
            while it != samples.end():
                dict[deref(it).first.c_str()] = string_map2dict(deref(it).second)
                inc(it)
            return dict

    property num_hom_ref:
        """ the count of homozygous for the ref allele gts"""
        def __get__(self):
            return self._thisptr.num_hom_ref

    property num_het:
        """ the count of heterozygous gts"""
        def __get__(self):
            return self._thisptr.num_het

    property num_hom_alt:
        """ the count of homozygous for the alt allele gts"""
        def __get__(self):
            return self._thisptr.num_hom_alt

    property num_unknown:
        """ the count of unknown gts"""
        def __get__(self):
            return self._thisptr.num_unknown

    property num_valid:
        """ the count of valid (i.e., !unknown) gts"""
        def __get__(self):
            return self._thisptr.num_valid

    property num_alt_alleles:
        """ the count of observed alternate alleles 
            (note alleles, not gts)
        """
        def __get__(self):
            return self._thisptr.num_alt_alleles

    property gts:
        """ return a list of the genotypes (AA, AG, GG) in order by sample"""
        def __get__(self):
            return string_vec2list(self._thisptr.gts)

    property gt_types:
        """ return a list of the genotypes (0, 1, 2) in order by sample"""
        def __get__(self):
            return int_vec2list(self._thisptr.gt_types)

    property gt_phases:
        """ return a list of the phase state of the genotypes 
            (True = phased, False = unphased) in order by sample
        """
        def __get__(self):
            return bool_vec2list(self._thisptr.gt_phases)

    property aaf:
        """return the alternate allele frequency for the variant"""
        def __get__(self):
            return self._thisptr.getAAF()

    property pi:
        """return the site-specific nucleotide diversity for the variant"""
        def __get__(self):
            return self._thisptr.getNucleotideDiversity()
