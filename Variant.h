#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <stdlib.h>
#include <assert.h>
#include <stack>
#include <queue>
#include "split.h"
#include "join.h"

using namespace std;

namespace vcf {

class Variant;

class VariantCallFile {

public:

    istream* file;

    string header;
    string line; // the current line
    string fileformat;
    string fileDate;
    string source;
    string reference;
    string phasing;
    map<string, string> infoTypes;
    map<string, int> infoCounts;
    map<string, string> formatTypes;
    map<string, int> formatCounts;
    vector<string> sampleNames;

    bool open(string& filename) {
        file = &_file;
        _file.open(filename.c_str(), ifstream::in);
        parsedHeader = parseHeader();
    }

    bool open(istream& stream) {
        file = &stream;
        parsedHeader = parseHeader();
    }

    bool open(ifstream& stream) {
        file = &stream;
        parsedHeader = parseHeader();
    }

    VariantCallFile(void) { }

    // open a file
    VariantCallFile(string& filename) : file(&_file) { 
        _file.open(filename.c_str(), ifstream::in);
        parsedHeader = parseHeader();
    }

    // use an existing stream as our file
    VariantCallFile(ifstream& stream) : file(&stream) { 
        parsedHeader = parseHeader();
    }

    VariantCallFile(istream& stream) : file(&stream) { 
        parsedHeader = parseHeader();
    }

    bool is_open(void) { return parsedHeader; }

    bool parseHeader(void);

    bool getNextVariant(Variant& var);

private:
    bool firstRecord;
    bool usingFile;
    ifstream _file;
    bool parsedHeader;

};

class Variant {

    friend ostream& operator<<(ostream& out, Variant& var);

public:

    VariantCallFile& vcf;
    string sequenceName;
    unsigned long position;
    string id;
    string ref;
    vector<string> alt;      // a list of all the alternate alleles present at this locus
    vector<string> alleles;  // a list all alleles (ref + alt) at this locus
                             // the indicies are organized such that the genotype codes (0,1,2,.etc.)
                             // correspond to the correct offest into the allelese vector.
                             // that is, alleles[0] = ref, alleles[1] = first alternate allele, etc.
    string filter;
    int quality;
    map<string, string> info;
    map<string, bool> infoFlags;
    vector<string> format;
    map<string, map<string, string> > samples;
    vector<string> sampleNames;

public:

    Variant(VariantCallFile& v)
        : sampleNames(v.sampleNames)
        , vcf(v)
    { }

    void parse(string& line);
    void addFilter(string& tag);
    bool getValueBool(string& key, string& sample);
    float getValueFloat(string& key, string& sample);
    string getValueString(string& key, string& sample);
    bool getSampleValueBool(string& key, string& sample);
    float getSampleValueFloat(string& key, string& sample);
    string getSampleValueString(string& key, string& sample);
    bool getInfoValueBool(string& key);
    float getInfoValueFloat(string& key);
    string getInfoValueString(string& key);

private:
    string lastFormat;

};

// from BamTools
// RuleToken implementation

struct RuleToken {

    // enums
    enum RuleTokenType { OPERAND = 0
                       , NUMBER
                       , BOOLEAN_VARIABLE
                       , NUMERIC_VARIABLE
                       , STRING_VARIABLE
                       , AND_OPERATOR
                       , OR_OPERATOR
                       , NOT_OPERATOR
                       , EQUAL_OPERATOR
                       , GREATER_THAN_OPERATOR
                       , LESS_THAN_OPERATOR
                       , LEFT_PARENTHESIS
                       , RIGHT_PARENTHESIS
                       };

    // constructor
    RuleToken(string token);
    RuleToken(void) 
        : type(BOOLEAN_VARIABLE)
        , state(false)
    { }

    // data members
    RuleTokenType type;
    string value;

    float number;
    string str;
    bool state;

    bool isVariable; // if this is a variable
    //bool isEvaluated; // when we evaluate variables

    RuleToken apply(RuleToken& other);

};

inline int priority(const RuleToken& token) {
    switch ( token.type ) {
        case ( RuleToken::NOT_OPERATOR )          : return 4;
        case ( RuleToken::EQUAL_OPERATOR )        : return 3;
        case ( RuleToken::GREATER_THAN_OPERATOR ) : return 3;
        case ( RuleToken::LESS_THAN_OPERATOR )    : return 3;
        case ( RuleToken::AND_OPERATOR )          : return 2;
        case ( RuleToken::OR_OPERATOR )           : return 1;
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
             token.type == RuleToken::LESS_THAN_OPERATOR
             );
}

inline bool isOperatorChar(const char& c) {
    return (c == '!' ||
            c == '&' ||
            c == '|' ||
            c == '=' ||
            c == '>' ||
            c == '<');
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

// converts the string into the specified type, setting r to the converted
// value and returning true/false on success or failure
template<typename T>
bool convert(const string& s, T& r) {
    istringstream iss(s);
    iss >> r;
    return (iss.fail() || iss.tellg() != s.size()) ? false : true;
}

void tokenizeFilterSpec(string& filterspec, stack<RuleToken>& tokens);


class VariantFilter {

public:

    enum VariantFilterType { SAMPLE = 0,
                             RECORD };

    string spec;
    queue<RuleToken> tokens; // tokens, infix notation
    queue<RuleToken> rules;  // tokens, prefix notation
    VariantFilterType type;
    VariantFilter(string filterspec, VariantFilterType filtertype);
    bool passes(Variant& var, string& sample);

};

} // end namespace VCF
