#ifndef BEDREADER_H
#define BEDREADER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <map>
#include <iterator>
#include <algorithm>
#include "intervaltree/IntervalTree.h"
#include "split.h"

using namespace std;

string strip(string const& str, char const* separators = " \t") {
    string::size_type const first = str.find_first_not_of(separators);
    return (first == string::npos) ? string()
        : str.substr(first, str.find_last_not_of(separators) - first + 1);
}

// stores the posiitional information of a bed target entry
class BedTarget {

public:

    string seq;  // sequence name
    int left;    // left position
    int right;   // right position, adjusted to 0-base
    string desc; // descriptive information, target name typically

    BedTarget(string s, int l, int r, string d = "")
        : seq(s)
        , left(l)
        , right(r)
        , desc(d)
    { }

};


class BedReader {

    bool _isOpen;
    ifstream file;

public:

    bool isOpen(void) { return _isOpen; }

    vector<BedTarget> targets;
    map<string, IntervalTree<BedTarget*> > intervals; // intervals by reference sequence

    vector<BedTarget> entries(void) {

        vector<BedTarget> entries;

        if (!isOpen()) {
            cerr << "bed targets file is not open" << endl;
            exit(1);
        }

        string line;
        while (std::getline(file, line)) {
            vector<string> fields = split(line, " \t");
            BedTarget entry(strip(fields[0]),
                            atoi(strip(fields[1]).c_str()),
                            atoi(strip(fields[2]).c_str()),
                            (fields.size() >= 4) ? strip(fields[3]) : "");
            entries.push_back(entry);
        }

        return entries;

    }

    vector<BedTarget*> targetsContained(BedTarget& target) {
        vector<Interval<BedTarget*> > results;
        intervals[target.seq].findContained(target.left, target.right, results);
        vector<BedTarget*> contained;
        for (vector<Interval<BedTarget*> >::iterator r = results.begin(); r != results.end(); ++r) {
            contained.push_back(r->value);
        }
        return contained;
    }

    vector<BedTarget*> targetsOverlapping(BedTarget& target) {
        vector<Interval<BedTarget*> > results;
        intervals[target.seq].findOverlapping(target.left, target.right, results);
        vector<BedTarget*> overlapping;
        for (vector<Interval<BedTarget*> >::iterator r = results.begin(); r != results.end(); ++r) {
            overlapping.push_back(r->value);
        }
        return overlapping;
    }

    BedReader(void)
	: _isOpen(false)
    { }

    BedReader(string& fname)
	: _isOpen(false) {
	open(fname);
    }

    void open(const string& fname) {
        file.open(fname.c_str());
	_isOpen = true;
        targets = entries();
        map<string, vector<Interval<BedTarget*> > > intervalsBySeq;
        for (vector<BedTarget>::iterator t = targets.begin(); t != targets.end(); ++t) {
            intervalsBySeq[t->seq].push_back(Interval<BedTarget*>(1 + t->left, t->right, &*t));
        }
        for (map<string, vector<Interval<BedTarget*> > >::iterator s = intervalsBySeq.begin(); s != intervalsBySeq.end(); ++s) {
            intervals[s->first] = IntervalTree<BedTarget*>(s->second);
        }
    }

};

#endif

