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
        

cdef extern from "Variant.h" namespace "vcf":
    cdef cppclass VariantCallFile:
        VariantCallFile()
        bint open(string &)
        bint getNextVariant(Variant)

    cdef cppclass Variant:
        Variant(VariantCallFile &)

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
        # custom methods
        float getAAF()
        float getNucleotideDiversity()
        
        bint parsedGenotypes
        void parseGenotypes()


# data structure conversion functions
cdef list string_vec2list(vector[string] sv):
    """
    convert an STL vector<string> to a list
    """
    cdef size_t size = sv.size(), i
    return [sv.at(i).c_str() for i in range(size)]

# data structure conversion functions
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

cdef class PyVariant:
    """
    A Python wrapper class for vcflib's Variant class.
    """
    cdef Variant *_thisptr
    def __init__(self):
        pass
    def __dealloc__(self):
        del self._thisptr

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
            if not self._thisptr.parsedGenotypes:
                self._thisptr.parseGenotypes()
            return self._thisptr.num_hom_ref
    property num_het:
        """ the count of heterozygous gts"""
        def __get__(self):
            if not self._thisptr.parsedGenotypes:
                self._thisptr.parseGenotypes()
            return self._thisptr.num_het
    property num_hom_alt:
        """ the count of homozygous for the alt allele gts"""
        def __get__(self):
            if not self._thisptr.parsedGenotypes:
                self._thisptr.parseGenotypes()
            return self._thisptr.num_hom_alt
    property num_unknown:
        """ the count of unknown gts"""
        def __get__(self):
            if not self._thisptr.parsedGenotypes:
                self._thisptr.parseGenotypes()
            return self._thisptr.num_unknown
    property gts:
        """ return a list of the genotypes (AA, AG) in order by sample"""
        def __get__(self):
            if not self._thisptr.parsedGenotypes:
                self._thisptr.parseGenotypes()
            return string_vec2list(self._thisptr.gts)
            
    # additional methods
    property aaf:
        """return the alternate allele frequency for the variant"""
        def __get__(self):
            return self._thisptr.getAAF()
    property pi:
        """return the site-specific nucleotide diversity for the variant"""
        def __get__(self):
            return self._thisptr.getNucleotideDiversity()



cdef PyVariant create_variant(Variant *v):
    cdef PyVariant pyvar = PyVariant.__new__(PyVariant)
    pyvar._thisptr = v
    return pyvar

cdef class VariantFile:
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
