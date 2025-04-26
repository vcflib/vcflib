/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"
#include <algorithm>
#include <getopt.h>
#include "convert.h"

using namespace std;
using namespace vcflib;

bool listContains(const list<string>& l, const string& v) {
    for (const auto& elm : l) {
        if (elm == v) return true;
    }
    return false;
}

void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] [vcf file]" << endl
         << endl
         << "Sorts the input (either stdin or file) using a streaming sort algorithm. Guarantees that the positional order is correct provided out-of-order variants are no more than 100 positions in the VCF file apart."
         << endl
         << "options:" << endl
         << endl
         << "    -h, --help             this dialog" << endl
         << "    -w, --window N         number of sites to sort (default 10000)" << endl
         << "    -a, --all              load all sites and then sort in memory" << endl;
    cerr << endl << "Type: transformation" << endl << endl;
}

int main(int argc, char** argv) {

    VariantCallFile variantFile;
    int sortSitesWindow = 10000;
    bool sortAll = false;

    int c;

    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
            {"window", required_argument, 0, 'w'},
            {"all", required_argument, 0, 'a'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "haw:",
                         long_options, &option_index);

        if (c == -1)
            break;

        string field;

        switch (c)
        {

        case 'w':
            if (!convert(optarg, sortSitesWindow)) {
                cerr << "could not parse --window, -w" << endl;
                exit(1);
            }
            break;

        case 'a':
            sortAll = true;
            break;

        case 'h':
            printSummary(argv);
            exit(0);
            break;

        default:
            break;
        }
    }

    if (optind == argc - 1) {
        string inputFilename = argv[optind];
        variantFile.open(inputFilename);
    } else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    cout << variantFile.header << endl;

    map<string, map<long int, map<string, vector<Variant> > > > records;
    long int back = 0;
    int numrecords = 0;
    list<string> sequenceNames;

    variantFile.parseSamples = false;
    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        //cerr << "at position " << var.sequenceName << ":" << var.position << endl;
        if (!listContains(sequenceNames, var.sequenceName)) {
            //cerr << "adding new sequence name " << var.sequenceName << endl;
            sequenceNames.push_back(var.sequenceName);
        }
        records[var.sequenceName][var.position][var.vrepr()].push_back(var);
        if (records[var.sequenceName][var.position].size() == 1) ++numrecords;
        if (!sortAll && numrecords > sortSitesWindow) {
            //cerr << "outputting a position" << endl;
            if (records[sequenceNames.front()].empty()) {
                //cerr << "end of reference sequence " << sequenceNames.front() << endl;
                sequenceNames.pop_front();
            }
            map<long int, map<string, vector<Variant> > >& frecords = records[sequenceNames.front()];
            map<string, vector<Variant> >& vars = frecords.begin()->second;
            for (map<string, vector<Variant> >::iterator v = vars.begin(); v != vars.end(); ++v) {
                for (vector<Variant>::iterator s = v->second.begin(); s != v->second.end(); ++s) {
                    cout << s->originalLine << endl;
                }
            }
            frecords.erase(frecords.begin());
            --numrecords;
        }
    }
    //cerr << "done processing input, cleaning up" << endl;
    for (const auto& sequenceName : sequenceNames) {
        map<long int, map<string, vector<Variant> > >& q = records[sequenceName];
        for (const auto& r : q) {
            for (const auto& v : r.second) {
                for (const auto& s : v.second) {
                    cout << s.originalLine << endl;
                }
            }
            --numrecords;
        }
    }
    //cerr << numrecords << " remain" << endl;

    return 0;

}
