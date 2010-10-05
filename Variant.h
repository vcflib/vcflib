#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <stdlib.h>
#include <assert.h>
#include "Split.h"

using namespace std;

class Variant;

class VariantCallFile : public ifstream {

public:

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

    bool openVCF(string& filename);

    bool getNextVariant(Variant& var);

private:
    bool firstRecord;

};

class Variant {

    friend ostream& operator<<(ostream& out, Variant& var);

public:

    VariantCallFile& vcf;
    string sequenceName;
    unsigned long position;
    string id;
    string ref;
    string alt;
    string filter;
    int quality;
    map<string, string> info;
    map<string, bool> infoFlags;
    vector<string> format;
    map<string, map<string, string> > samples;
    vector<string> sampleNames;

public:

    Variant(vector<string>& sns, VariantCallFile& v)
        : sampleNames(sns)
        , vcf(v)
    { }

    void parse(string& line);

private:
    string lastFormat;

};

// op type
enum VariantFilterOpType {
    FILTER_GREATER_THAN_OR_EQUAL,
    FILTER_LESS_THAN_OR_EQUAL,
    FILTER_EQUAL,
    FILTER_NOT_EQUAL,
    FILTER_GREATER_THAN,
    FILTER_LESS_THAN,
    FILTER_FLAG
};

template <class T>
bool applyFilter(VariantFilterOpType op, T data, T value) {
    switch (op) {
        case FILTER_GREATER_THAN_OR_EQUAL:
            return data >= value;
            break;
        case FILTER_LESS_THAN_OR_EQUAL:
            return data <= value;
            break;
        case FILTER_GREATER_THAN:
            return data > value;
            break;
        case FILTER_LESS_THAN:
            return data < value;
            break;
        case FILTER_EQUAL:
            return data == value;
            break;
        case FILTER_NOT_EQUAL:
            return data != value;
            break;
        default:
            break;
    }
    return false;
}

class VariantFilter {

public:

    string spec;
    string field;   // field ID in the variant record
    bool isInfoField;  // true if an info field, 
    VariantFilterOpType op;
    string value;  // value against which to check the key

    VariantFilter(string filterspec, bool infof = true);
    bool passes(Variant& var);

};

