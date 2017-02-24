#ifndef __VARIANT_H
#define __VARIANT_H

#include <vector>
#include <list>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <stdlib.h>
#include <assert.h>
#include <stack>
#include <queue>
#include <set>
#include <functional>
#include <cstdio>
#include "split.h"
#include "join.h"
#include "tabix.hpp"
#include "SmithWatermanGotoh.h"
#include "disorder.h"
#include "ssw_cpp.hpp"
#include "convert.h"
#include "multichoose.h"
#include "Fasta.h"
extern "C" {
    #include "filevercmp.h"
}

using namespace std;

namespace vcflib {

class Variant;

enum VariantFieldType { FIELD_FLOAT = 0
                      , FIELD_INTEGER
                      , FIELD_BOOL
                      , FIELD_STRING
                      , FIELD_UNKNOWN
                      };

enum VariantFieldNumber { ALLELE_NUMBER = -2
                        , GENOTYPE_NUMBER = -1
                        };

const int INDEX_NONE = -1;
const int NULL_ALLELE = -1;

VariantFieldType typeStrToFieldType(string& typeStr);
ostream& operator<<(ostream& out, VariantFieldType type);

typedef map<string, map<string, vector<string> > > Samples;
typedef vector<pair<int, string> > Cigar;

class VariantCallFile {

public:

    istream* file;
    Tabix* tabixFile;

    bool usingTabix;
    string vcf_header;


    string header;
    string line; // the current line
    string fileformat;
    string fileDate;
    string source;
    string reference;
    string phasing;
    map<string, VariantFieldType> infoTypes;
    map<string, int> infoCounts;
    map<string, VariantFieldType> formatTypes;
    map<string, int> formatCounts;
    vector<string> sampleNames;
    bool parseSamples;
    bool _done;

    void updateSamples(vector<string>& newSampleNames);
    string headerWithSampleNames(vector<string>& newSamples); // non-destructive, for output
    void addHeaderLine(string line);
    void removeInfoHeaderLine(string line);
    void removeGenoHeaderLine(string line);
    vector<string> infoIds(void);
    vector<string> formatIds(void);

    bool open(string& filename) {
        vector<string> filenameParts = split(filename, ".");
        if (filenameParts.back() == "gz" || filenameParts.back() == "bgz") {
            return openTabix(filename);
        } else {
            return openFile(filename);
        }
    }

    bool openFile(string& filename) {
        file = &_file;
        _file.open(filename.c_str(), ifstream::in);
        parsedHeader = parseHeader();
        return parsedHeader;
    }

    bool openTabix(string& filename) {
        usingTabix = true;
        tabixFile = new Tabix(filename);
        parsedHeader = parseHeader();
        return parsedHeader;
    }

    bool open(istream& stream) {
        file = &stream;
        parsedHeader = parseHeader();
        return parsedHeader;
    }

    bool open(ifstream& stream) {
        file = &stream;
        parsedHeader = parseHeader();
        return parsedHeader;
    }

    bool openForOutput(string& headerStr) {
        parsedHeader = parseHeader(headerStr);
        return parsedHeader;
    }

VariantCallFile(void) : usingTabix(false), parseSamples(true), justSetRegion(false), parsedHeader(false) { }
    ~VariantCallFile(void) {
        if (usingTabix) {
            delete tabixFile;
        }
    }

    bool is_open(void) { return parsedHeader; }

    bool eof(void) { return _file.eof(); }

    bool done(void) { return _done; }

    bool parseHeader(string& headerStr);

    bool parseHeader(void);

    bool getNextVariant(Variant& var);

    bool setRegion(string region);
    bool setRegion(string seq, long int start, long int end = 0);
    vector<string> getHeaderLinesFromFile();

private:
    bool firstRecord;
    bool justSetRegion;
    bool usingFile;
    ifstream _file;
    bool parsedHeader;

};

class VariantAllele {
    friend ostream& operator<<(ostream& out, VariantAllele& var);
    friend bool operator<(const VariantAllele& a, const VariantAllele& b);
    friend VariantAllele operator+(const VariantAllele& a, const VariantAllele& b);
public:
    string ref;
    string alt;
    string repr;
    long position;
    /* // TODO
    bool isSNP(void);
    bool isMNP(void);
    bool isInsertion(void);
    bool isDeletion(void);
    bool isIndel(void);
    */
    VariantAllele(string r, string a, long p)
        : ref(r), alt(a), position(p)
    {
        stringstream s;
        s << position << ":" << ref << "/" << alt;
        repr = s.str();
    }
};

class Variant {

    friend ostream& operator<<(ostream& out, Variant& var);
    
public:

    string sequenceName;
    long position;
    long zeroBasedPosition(void);
    string id;
    string ref;
    vector<string> alt;      // a list of all the alternate alleles present at this locus
    vector<string> alleles;  // a list all alleles (ref + alt) at this locus
                             // the indicies are organized such that the genotype codes (0,1,2,.etc.)
                             // correspond to the correct offest into the allelese vector.
                             // that is, alleles[0] = ref, alleles[1] = first alternate allele, etc.
    string vrepr(void);  // a comparable record of the variantion described by the record
    set<string> altSet(void);  // set of alleles, rather than vector of them
    map<string, int> altAlleleIndexes;  // reverse lookup for alleles
    map<string, vector<VariantAllele> > parsedAlternates(bool includePreviousBaseForIndels = false,
                                                         bool useMNPs = false,
                                                         bool useEntropy = false,
                                                         float matchScore = 10.0f,
                                                         float mismatchScore = -9.0f,
                                                         float gapOpenPenalty = 15.0f,
                                                         float gapExtendPenalty = 6.66f,
                                                         float repeatGapExtendPenalty = 0.0f,
                                                         string flankingRefLeft = "",
                                                         string flankingRefRight = "");
    // the same output format as parsedAlternates, without parsing
    map<string, vector<VariantAllele> > flatAlternates(void);

    map<string, string> extendedAlternates(long int newPosition, long int length);

    // Convert a structural variant the canonical VCF4.2 format using a reference.
    // returns true if the variant is canonicalized, false otherwise.
    bool canonicalize_sv(FastaReference& ref, vector<FastaReference*> insertions, int interval_sz = -1);

    bool is_sv();

    string originalLine; // the literal of the record, as read
    // TODO
    // the ordering of genotypes for the likelihoods is given by: F(j/k) = (k*(k+1)/2)+j
    // vector<pair<int, int> > genotypes;  // indexes into the alleles, ordered as per the spec
    string filter;
    double quality;
    VariantFieldType infoType(string& key);
    map<string, vector<string> > info;  // vector<string> allows for lists by Genotypes or Alternates
    map<string, bool> infoFlags;
    VariantFieldType formatType(string& key);
    vector<string> format;
    map<string, map<string, vector<string> > > samples;  // vector<string> allows for lists by Genotypes or Alternates
    vector<string> sampleNames;
    vector<string> outputSampleNames;
    VariantCallFile* vcf;

    //void addInfoInt(string& tag, int value);
    //void addInfoFloat(string& tag, double value);
    //void addInfoString(string& tag, string& value);

    void removeAlt(string& altallele);

public:

    Variant() { }

    Variant(VariantCallFile& v)
        : sampleNames(v.sampleNames)
        , outputSampleNames(v.sampleNames)
        , vcf(&v)
    { }

    void setVariantCallFile(VariantCallFile& v);
    void setVariantCallFile(VariantCallFile* v);

    void parse(string& line, bool parseSamples = true);
    void addFilter(string& tag);
    bool getValueBool(string& key, string& sample, int index = INDEX_NONE);
    double getValueFloat(string& key, string& sample, int index = INDEX_NONE);
    string getValueString(string& key, string& sample, int index = INDEX_NONE);
    bool getSampleValueBool(string& key, string& sample, int index = INDEX_NONE);
    double getSampleValueFloat(string& key, string& sample, int index = INDEX_NONE);
    string getSampleValueString(string& key, string& sample, int index = INDEX_NONE);
    bool getInfoValueBool(string& key, int index = INDEX_NONE);
    double getInfoValueFloat(string& key, int index = INDEX_NONE);
    string getInfoValueString(string& key, int index = INDEX_NONE);
    void printAlt(ostream& out);      // print a comma-sep list of alternate alleles to an ostream
    void printAlleles(ostream& out);  // print a comma-sep list of *all* alleles to an ostream
    int getAltAlleleIndex(string& allele);
    void updateAlleleIndexes(void);
    void addFormatField(string& key);
    void setOutputSampleNames(vector<string>& outputSamples);
    map<pair<int, int>, int> getGenotypeIndexesDiploid(void);
    int getNumSamples(void);
    int getNumValidGenotypes(void);
    string getGenotype(string& sample);
    bool isPhased(void);
    // TODO
    //void setInfoField(string& key, string& val);

private:

    string lastFormat;

};


// from BamTools
// RuleToken implementation

class RuleToken {

public:

    // enums
    enum RuleTokenType { OPERAND = 0
                       , NUMBER
                       , BOOLEAN_VARIABLE
                       , NUMERIC_VARIABLE
                       , STRING_VARIABLE
                       , AND_OPERATOR
                       , OR_OPERATOR
                       , ADD_OPERATOR
                       , SUBTRACT_OPERATOR
                       , MULTIPLY_OPERATOR
                       , DIVIDE_OPERATOR
                       , NOT_OPERATOR
                       , EQUAL_OPERATOR
                       , GREATER_THAN_OPERATOR
                       , LESS_THAN_OPERATOR
                       , LEFT_PARENTHESIS
                       , RIGHT_PARENTHESIS
                       };

    // constructor
    RuleToken(string token, map<string, VariantFieldType>& variables);
    RuleToken(void) 
        : type(BOOLEAN_VARIABLE)
        , state(false)
    { }

    // data members
    RuleTokenType type;
    string value;

    double number;
    string str;
    bool state;

    bool isVariable; // if this is a variable
    //bool isEvaluated; // when we evaluate variables

    RuleToken apply(RuleToken& other);

};

inline int priority(const RuleToken& token) {
    switch ( token.type ) {
        case ( RuleToken::MULTIPLY_OPERATOR )     : return 8;
        case ( RuleToken::DIVIDE_OPERATOR )       : return 8;
        case ( RuleToken::ADD_OPERATOR )          : return 7;
        case ( RuleToken::SUBTRACT_OPERATOR )     : return 7;
        case ( RuleToken::NOT_OPERATOR )          : return 6;
        case ( RuleToken::EQUAL_OPERATOR )        : return 5;
        case ( RuleToken::GREATER_THAN_OPERATOR ) : return 5;
        case ( RuleToken::LESS_THAN_OPERATOR )    : return 5;
        case ( RuleToken::AND_OPERATOR )          : return 4;
        case ( RuleToken::OR_OPERATOR )           : return 3;
        case ( RuleToken::LEFT_PARENTHESIS )      : return 0;
        case ( RuleToken::RIGHT_PARENTHESIS )     : return 0;
        default: cerr << "invalid token type" << endl; exit(1);
    }
}

inline bool isRightAssociative(const RuleToken& token) {
    return (token.type == RuleToken::NOT_OPERATOR ||
            token.type == RuleToken::LEFT_PARENTHESIS);
}

inline bool isLeftAssociative(const RuleToken& token) {
    return !isRightAssociative(token);
}

inline bool isLeftParenthesis(const RuleToken& token) {
    return ( token.type == RuleToken::LEFT_PARENTHESIS );
}

inline bool isRightParenthesis(const RuleToken& token) {
    return ( token.type == RuleToken::RIGHT_PARENTHESIS );
}

inline bool isOperand(const RuleToken& token) {
    return ( token.type == RuleToken::OPERAND || 
             token.type == RuleToken::NUMBER ||
             token.type == RuleToken::NUMERIC_VARIABLE ||
             token.type == RuleToken::STRING_VARIABLE ||
             token.type == RuleToken::BOOLEAN_VARIABLE
           );
}

inline bool isOperator(const RuleToken& token) {
    return ( token.type == RuleToken::AND_OPERATOR ||
             token.type == RuleToken::OR_OPERATOR  ||
             token.type == RuleToken::NOT_OPERATOR ||
             token.type == RuleToken::EQUAL_OPERATOR ||
             token.type == RuleToken::GREATER_THAN_OPERATOR ||
             token.type == RuleToken::LESS_THAN_OPERATOR ||
             token.type == RuleToken::MULTIPLY_OPERATOR ||
             token.type == RuleToken::DIVIDE_OPERATOR ||
             token.type == RuleToken::ADD_OPERATOR ||
             token.type == RuleToken::SUBTRACT_OPERATOR
             );
}

inline bool isOperatorChar(const char& c) {
    return (c == '!' ||
            c == '&' ||
            c == '|' ||
            c == '=' ||
            c == '>' ||
            c == '<' ||
            c == '*' ||
            c == '/' ||
            c == '+' ||
            c == '-');
}

inline bool isParanChar(const char& c) {
    return (c == '(' || c == ')');
}

inline bool isNumeric(const RuleToken& token) {
    return token.type == RuleToken::NUMERIC_VARIABLE;
}

inline bool isString(const RuleToken& token) {
    return token.type == RuleToken::STRING_VARIABLE;
}

inline bool isBoolean(const RuleToken& token) {
    return token.type == RuleToken::BOOLEAN_VARIABLE;
}

inline bool isVariable(const RuleToken& token) {
    return isNumeric(token) || isString(token) || isBoolean(token);
}

void tokenizeFilterSpec(string& filterspec, stack<RuleToken>& tokens, map<string, VariantFieldType>& variables);


class VariantFilter {

public:

    enum VariantFilterType { SAMPLE = 0,
                             RECORD };

    string spec;
    queue<RuleToken> tokens; // tokens, infix notation
    queue<RuleToken> rules;  // tokens, prefix notation
    VariantFilterType type;
    VariantFilter(string filterspec, VariantFilterType filtertype, map<string, VariantFieldType>& variables);
    bool passes(Variant& var, string& sample); // all alts pass
    bool passes(Variant& var, string& sample, string& allele);
    void removeFilteredGenotypes(Variant& var, bool keepInfo);

};


// genotype manipulation

// TODO
//map<string, int> decomposeGenotype(string& genotype);

vector<int> decomposePhasedGenotype(const string& genotype);
map<int, int> decomposeGenotype(const string& genotype);

string genotypeToString(const map<int, int>& genotype);

string phasedGenotypeToString(const vector<int>& genotype);

bool isHet(const map<int, int>& genotype);

bool isHom(const map<int, int>& genotype);

bool hasNonRef(const map<int, int>& genotype);

bool isHomRef(const map<int, int>& genotype);

bool isHomNonRef(const map<int, int>& genotype);

bool isNull(const map<int, int>& genotype);

int ploidy(const map<int, int>& genotype);

string unionInfoHeaderLines(string& s1, string& s2);

// genotype likelihood ordering

list<list<int> > glorder(int ploidy, int alleles);
list<list<int> > _glorder(int ploidy, int alleles);
list<int> glsWithAlt(int alt, int ploidy, int numalts);
map<int, int> glReorder(int ploidy, int numalts, map<int, int>& alleleIndexMapping, vector<int>& altsToRemove);

vector<string>& unique(vector<string>& strings);

string varCigar(vector<VariantAllele>& vav, bool xForMismatch = false);
string mergeCigar(const string& c1, const string& c2);
vector<pair<int, string> > splitCigar(const string& cigarStr);
list<pair<int, string> > splitCigarList(const string& cigarStr);
int cigarRefLen(const vector<pair<int, char> >& cigar);
int cigarRefLen(const vector<pair<int, string> >& cigar);
vector<pair<int, string> > cleanCigar(const vector<pair<int, string> >& cigar);
string joinCigar(const vector<pair<int, string> >& cigar);
string joinCigar(const vector<pair<int, char> >& cigar);
string joinCigarList(const list<pair<int, string> >& cigar);
bool isEmptyCigarElement(const pair<int, string>& elem);

// for sorting, generating maps ordered by chromosome name
class ChromNameCompare {
public:
    bool operator()(const string& a, const string& b) const {
        return (filevercmp(a.c_str(), b.c_str()) < 0);
    }
};

class VCFHeader
{
public:
    VCFHeader();
    ~VCFHeader() {}

    /*
     * Adds header_column to this->header_columns if
     * it doesn't already exits.
     */
    void addHeaderColumn(const string& header_column);

    /*
     * Adds meta_line to either header_lines or header_lists.
     *
     * We parse out the ##_type_ from meta_line
     * - If the meta_line ##_type_ is a key in header_lines then meta_line is added to header_lines
     * - If the meta_line ##_type_ is a key in header_lists then meta_line is added to header_lists[##_type_] vector<string>
     *    Unless that header_lists[##_type_] vector already contains the ID that is in meta_line, in that case it is not added
     */
    void addMetaInformationLine(const string& meta_line);

    /*
     * Converts header_lines, header_lists and header_columns to a proper VCF header
     */
    string getHeaderString();

private:
    VCFHeader(const VCFHeader& vcfHeader); // Do not implement the copy constructor, there is no reason to add this functionality
    VCFHeader& operator=(const VCFHeader& vcfHeader); // Do not implement operator=, there is no reason to add this functionality

    /*
     * This is a helper function that determines if the ID substring contained in meta_line
     * exists as a ID substring within the vector<string> meta_lines. Returns true if
     * the ID exists within the vector and false otherwise.
     */
    bool metaInfoIdExistsInVector(const string& meta_line, vector<string>& meta_lines);

    /*
     * header_line_names_ordered contains all the header lines that
     * are available and in the expected order for a valid VCF file
     */
    vector<string> header_line_names_ordered;
    /*
     * header_list_names_ordered contains all the header lists that
     * are available and in the expected order for a valid VCF file
     */
    vector<string> header_list_names_ordered;

    /*
     * header_columns is set by the constructor to contain the 8 manditory VCF fields.
     * Also, unique header_columns for each of the vcf files are added as well.
     * Duplicates are not allowed, to prevent duplicates use addHeaderColumn when adding header columns
     */
    vector<string> header_columns;

    /* 
     * the maps we're going to be using will be case-insensitive
     * so that "fileFormat" and "fileformat" hash to the same item.
     */
    struct stringcasecmp : binary_function<string, string, bool> {
        struct charcasecmp : public std::binary_function<unsigned char, unsigned char, bool> {
            bool operator() (const unsigned char& c1, const unsigned char& c2) const {
                return tolower (c1) < tolower (c2); 
            }
        };
        bool operator() (const std::string & s1, const std::string & s2) const {
            return std::lexicographical_compare(s1.begin(), s1.end(), s2.begin(), s2.end(), charcasecmp());
        }
    };

    // contains all the ##_types_ as keys, the value is either empty or a VCF file has set it
    map<string, string, stringcasecmp> header_lines; 

    // contains all the ##_types_ as keys, the value is a vector of ##_type_ (since there can be duplicate #INFO for example, duplicate ids are not allowed)
    map<string, vector<string>, stringcasecmp> header_lists; 

};

} // end namespace VCF

#endif
