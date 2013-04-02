#include "Variant.h"
#include <algorithm>

using namespace std;
using namespace vcf;

bool listContains(list<string>& l, string& v) {
    for (list<string>::iterator i = l.begin(); i != l.end(); ++i) {
        if (*i == v) return true;
    }
    return false;
}

int main(int argc, char** argv) {

    VariantCallFile variantFile;

    if (argc > 1) {
        string filename = argv[1];
        variantFile.open(filename);
    } else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    cout << variantFile.header << endl;

    map<string, map<long int, vector<Variant> > > records;
    long int back = 0;
    int sortSitesWindow = 100;
    int numrecords = 0;
    list<string> sequenceNames;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        //cerr << "at position " << var.sequenceName << ":" << var.position << endl;
        if (!listContains(sequenceNames, var.sequenceName)) {
            //cerr << "adding new sequence name " << var.sequenceName << endl;
            sequenceNames.push_back(var.sequenceName);
        }
        records[var.sequenceName][var.position].push_back(var);
        if (records[var.sequenceName][var.position].size() == 1) ++numrecords;
        if (numrecords > sortSitesWindow) {
            //cerr << "outputting a position" << endl;
            if (records[sequenceNames.front()].empty()) {
                //cerr << "end of reference sequence " << sequenceNames.front() << endl;
                sequenceNames.pop_front();
            }
            map<long int, vector<Variant> >& frecords = records[sequenceNames.front()];
            vector<Variant>& vars = frecords.begin()->second;
            for (vector<Variant>::iterator v = vars.begin(); v != vars.end(); ++v) {
                cout << *v << endl;
            }
            frecords.erase(frecords.begin());
            --numrecords;
        }
    }
    //cerr << "done processing input, cleaning up" << endl;
    for (list<string>::iterator s = sequenceNames.begin(); s != sequenceNames.end(); ++s) {
        map<long int, vector<Variant> >& q = records[*s];
        for (map<long int, vector<Variant> >::iterator r = q.begin(); r != q.end(); ++r) {
            for (vector<Variant>::iterator v = r->second.begin(); v != r->second.end(); ++v) {
                cout << *v << endl;
            }
            --numrecords;
        }
    }
    //cerr << numrecords << " remain" << endl;

    return 0;

}

