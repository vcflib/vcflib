from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.pair cimport pair


cdef extern from "split.h":
    # split a string on a single delimiter character (delim)
    vector[string]& split(const string &s, char delim, vector[string] &elems)
    vector[string]  split(const string &s, char delim)
    # split a string on any character found in the string of delimiters (delims)
    vector[string]& split(const string &s, const string& delims, vector[string] &elems)
    vector[string]  split(const string &s, const string& delims)


cdef extern from "Variant.h" namespace "vcf":

    cdef enum VariantFieldType:
        FIELD_FLOAT = 0
        FIELD_INTEGER
        FIELD_BOOL
        FIELD_STRING
        FIELD_UNKNOWN
    
    cdef enum VariantFieldNumber:
        ALLELE_NUMBER = -2
        GENOTYPE_NUMBER = -1
    
    const int INDEX_NONE = -1
    const int NULL_ALLELE = -1
    
    VariantFieldType typeStrToFieldType(string& typeStr)

    cdef cppclass VariantCallFile:
#        istream* file
#        Tabix* tabixFile
        bool usingTabix
        string header
        string line
        string fileformat
        string fileDate
        string source
        string reference
        string phasing
        map[string, VariantFieldType] infoTypes
        map[string, int] infoCounts
        map[string, VariantFieldType] formatTypes
        map[string, int] formatCounts
        vector[string] sampleNames
        bool parseSamples
        bool _done
        void updateSamples(vector[string]& newSampleNames)
        void addHeaderLine(string line)
        void removeInfoHeaderLine(string line)
        void removeGenoHeaderLine(string line)
        vector[string] infoIds()
        vector[string] formatIds()
        vector[string] filterIds()
        bool open(string& filename) 
        bool openFile(string& filename)
        bool openTabix(string& filename)
#        bool open(istream& stream) 
#        bool open(ifstream& stream) 
        bool openForOutput(string& headerStr)
        bool is_open()
        bool eof()
        bool done()
        bool parseHeader(string& headerStr)
        bool parseHeader()
        bool getNextVariant(Variant& var)
        bool setRegion(string region)
        bool setRegion(string seq, long int start, long int end)


    cdef cppclass VariantAllele:
        string ref
        string alt
        string repr
        long position
        VariantAllele(string r, string a, long p)
        
    
    cdef cppclass Variant:
        string sequenceName
        long position
        string id
        string ref
        vector[string] alt
        vector[string] alleles
        map[string, int] altAlleleIndexes
        map[string, vector[VariantAllele] ] parsedAlternates(bool includePreviousBaseForIndels,
                                 bool useMNPs,
                                 bool useEntropy,
                                 float matchScore,
                                 float mismatchScore,
                                 float gapOpenPenalty,
                                 float gapExtendPenalty,
                                 float repeatGapExtendPenalty,
                                 string flankingRefLeft,
                                 string flankingRefRight)
        map[string, string] extendedAlternates(long int newPosition, long int length)
        string originalLine
        string filter
        double quality
        VariantFieldType infoType(string& key)
        map[string, vector[string]] info
        map[string, bool] infoFlags
        VariantFieldType formatType(string& key)
        vector[string] format
        map[string, map[string, vector[string]]] samples
        vector[string] sampleNames
        vector[string] outputSampleNames
        VariantCallFile* vcf
        void removeAlt(string& altallele)
        Variant()
        Variant(VariantCallFile& v)
        void setVariantCallFile(VariantCallFile& v)
        void setVariantCallFile(VariantCallFile* v)
        void parse(string& line, bool parseSamples)
        void addFilter(string& tag)
        bool getValueBool(string& key, string& sample, int index)
        double getValueFloat(string& key, string& sample, int index)
        string getValueString(string& key, string& sample, int index)
        bool getSampleValueBool(string& key, string& sample, int index)
        double getSampleValueFloat(string& key, string& sample, int index)
        string getSampleValueString(string& key, string& sample, int index)
        bool getInfoValueBool(string& key, int index)
        double getInfoValueFloat(string& key, int index)
        string getInfoValueString(string& key, int index)
#        void printAlt(ostream& out)
#        void printAlleles(ostream& out)
        int getAltAlleleIndex(string& allele)
        void updateAlleleIndexes()
        void addFormatField(string& key)
        void setOutputSampleNames(vector[string]& outputSamples)
        map[pair[int, int], int] getGenotypeIndexesDiploid()
        int getNumSamples()
        int getNumValidGenotypes()
        

