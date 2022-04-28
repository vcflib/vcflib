#pragma once

#include "convert.h"
#include "join.h"
#include "split.h"
#include "allele.hpp"
#include <vector>
#include <set>
#include <iostream>

namespace vcflib {

using namespace std;

typedef vector<pair<int, string> > Cigar;
string varCigar(vector<VariantAllele>& vav, bool xForMismatch = false);
string mergeCigar(const string& c1, const string& c2);
vector<pair<int, string> > splitUnpackedCigar(const string& cigarStr);
vector<pair<int, string> > splitCigar(const string& cigarStr);
list<pair<int, string> > splitCigarList(const string& cigarStr);
int cigarRefLen(const vector<pair<int, char> >& cigar);
int cigarRefLen(const vector<pair<int, string> >& cigar);
vector<pair<int, string> > cleanCigar(const vector<pair<int, string> >& cigar);
string joinCigar(const vector<pair<int, string> >& cigar);
string joinCigar(const vector<pair<int, char> >& cigar);
string joinCigarList(const list<pair<int, string> >& cigar);
bool isEmptyCigarElement(const pair<int, string>& elem);

}
