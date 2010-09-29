#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <stdlib.h>
#include "Split.h"

using namespace std;


class Variant {

public:

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

    Variant(vector<string>& sns)
        : sampleNames(sns)
    { }


    void parse(string& line) {

        // clean up potentially variable data structures
        info.clear();
        format.clear();

        // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT [SAMPLE1 .. SAMPLEN]
        vector<string> fields = split(line, '\t');
        sequenceName = fields.at(0);
        char* end; // dummy variable for strtoll
        position = strtoll(fields.at(1).c_str(), &end, 10);
        id = fields.at(2);
        ref = fields.at(3);
        alt = fields.at(4); // TODO handle multi-allelic situations
        quality = atoi(fields.at(5).c_str());
        filter = fields.at(6);
        vector<string> infofields = split(fields.at(7), ';');
        for (vector<string>::iterator f = infofields.begin(); f != infofields.end(); ++f) {
            vector<string> kv = split(*f, '=');
            if (kv.size() == 2) {
                info[kv.at(0)] = kv.at(1);
            } else if (kv.size() == 1) {
                infoFlags[kv.at(0)] = true;
            }
        }
        format = split(fields.at(8), ':');
        if (fields.at(8) != lastFormat) {
            samples.clear();
            lastFormat = fields.at(8);
        }
        vector<string>::iterator sampleName = sampleNames.begin();
        for (vector<string>::iterator sample = fields.begin() + 9; sample != fields.end(); ++sample, ++sampleName) {
            if (*sample == ".") {
                samples[*sampleName].clear();
                continue;
            }
            vector<string> infofields = split(*sample, ':');
            vector<string>::iterator i = infofields.begin();
            for (vector<string>::iterator f = format.begin(); f != format.end(); ++f) {
                samples[*sampleName][*f] = *i; ++i;
            }
        }
        //return true; // we should be catching exceptions...
    }

    friend ostream& operator<<(ostream& out, Variant& var);

private:
    string lastFormat; // last format string, so we don't have to rebuild the sample map on every step

};

ostream& operator<<(ostream& out, Variant& var) {
    out << var.sequenceName << "\t"
        << var.position << "\t"
        << var.id << "\t"
        << var.ref << "\t"
        << var.alt << "\t"
        << var.filter << "\t"
        << var.quality << "\t";
    for (map<string, string>::iterator i = var.info.begin(); i != var.info.end(); ++i) {
        out << ((i == var.info.begin()) ? "" : ";") << i->first << "=" << i->second;
    }
    for (map<string, bool>::iterator i = var.infoFlags.begin(); i != var.infoFlags.end(); ++i) {
        out << ((i == var.infoFlags.end()) ? "" : ";") << i->first;
    }
    out << "\t";
    for (vector<string>::iterator f = var.format.begin(); f != var.format.end(); ++f) {
        out << ((f == var.format.begin()) ? "" : ":") << *f;
    }
    for (vector<string>::iterator s = var.sampleNames.begin(); s != var.sampleNames.end(); ++s) {
        out << "\t";
        map<string, map<string, string> >::iterator sampleItr = var.samples.find(*s);
        if (sampleItr == var.samples.end()) {
            out << ".";
        } else {
            map<string, string>& sample = sampleItr->second;
            if (sample.size() == 0) {
                out << ".";
            } else {
                for (vector<string>::iterator f = var.format.begin(); f != var.format.end(); ++f) {
                    out << ((f == var.format.begin()) ? "" : ":") << sample[*f];
                }
            }
        }
    }
    return out;
}


class VariantCallFile : public ifstream {

public:

    string header;
    string line; // the current line
    string fileformat;
    string fileDate;
    string source;
    string reference;
    string phasing;
    vector<string> infoDescriptions;
    map<string, pair<int, string> > infoTypes;
    vector<string> formatDescriptions;
    map<string, pair<int, string> > formatTypes;
    vector<string> sampleNames;

    bool openVCF(string& filename) {
        open(filename.c_str(), ifstream::in);
        if (!is_open()) {
            cerr << "could not open " << filename << endl;
            return false;
        }
        header = "";
        while (std::getline(*this, line)) {
            if (line.substr(0,2) == "##") {
                // meta-information lines
                // TODO parse
            } else if (line.substr(0,1) == "#") {
                // field name line
                vector<string> fields = split(line, '\t');
                sampleNames.resize(fields.size() - 9);
                copy(fields.begin() + 9, fields.end(), sampleNames.begin());
            } else {
                // done with header
                firstRecord = true;
                header.resize(header.size() - 1);
                return true;
            }
            header += line + '\n';
        }
    }

    bool getNextVariant(Variant& var) {
        if (firstRecord) {
            var.parse(line);
            firstRecord = false;
            return true;
        }
        if (std::getline(*this, line)) {
            var.parse(line);
            return true;
        } else {
            return false;
        }
    }

private:
    bool firstRecord;

};
