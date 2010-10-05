#include "Variant.h"


void Variant::parse(string& line) {

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

ostream& operator<<(ostream& out, Variant& var) {
    out << var.sequenceName << "\t"
        << var.position << "\t"
        << var.id << "\t"
        << var.ref << "\t"
        << var.alt << "\t"
        << var.quality << "\t"
        << var.filter << "\t";
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


VariantFilter::VariantFilter(string filterspec, bool infof) {
    isInfoField = infof; // marks if this filters on info-items or all genotypes
    size_t foundKey = filterspec.find(' ');
    if (foundKey == string::npos) {
        cerr << "could not parse filter spec " << filterspec << endl;
        exit(1);
    } else {
        field = filterspec.substr(0, foundKey);
        cerr << "found field " << field << endl;
        size_t filterOp;
        if ((filterOp = filterspec.find("<=")) != string::npos) {
            cerr << "found <=" << endl;
            op = FILTER_LESS_THAN_OR_EQUAL;
        } else if ((filterOp = filterspec.find(">=")) != string::npos) {
            cerr << "found >=" << endl;
            op = FILTER_GREATER_THAN_OR_EQUAL;
        } else if ((filterOp = filterspec.find("==")) != string::npos) {
            cerr << "found ==" << endl;
            op = FILTER_EQUAL;
        } else if ((filterOp = filterspec.find("!=")) != string::npos) {
            cerr << "found !=" << endl;
            op = FILTER_NOT_EQUAL;
        } else if ((filterOp = filterspec.find("<")) != string::npos) {
            cerr << "found <" << endl;
            op = FILTER_LESS_THAN;
        } else if ((filterOp = filterspec.find(">")) != string::npos) {
            cerr << "found >" << endl;
            op = FILTER_GREATER_THAN;
        } else {
            cerr << "regarding as boolean flag filter" << endl;
            op = FILTER_FLAG;
        }
        value = filterspec.substr(filterOp);
        cerr << "got value " << value << endl;
    }
}

bool VariantFilter::passes(Variant& var) {
    if (isInfoField) {
        string type = var.vcf.infoTypes[field];
        if (type == "Flag") {
            map<string, bool>::iterator f = var.infoFlags.find(field);
            if (f != var.infoFlags.end()) {
                return true;
            } else {
                return false;
            }
        } else if (type == "Integer") {
            return applyFilter(op, atoi(var.info[field].c_str()), atoi(value.c_str()));
        } else if (type == "Float") {
            return applyFilter(op, atof(var.info[field].c_str()), atof(value.c_str()));
        } else if (type == "String") {
            return applyFilter(op, var.info[field], value);
        } else {
            cerr << "unsupported field type: " << type << endl;
            exit(1);
        }
    } else { // go through all samples
        string type = var.vcf.formatTypes[field];
        for (vector<string>::iterator s = var.sampleNames.begin(); s != var.sampleNames.end(); ++s) {
            map<string, string>& sample = var.samples[*s];
            // unclear what a format-field flag would be
            /*if (type == "Flag") {
                map<string, string>::iterator f = var.infoFlags.find(key);
                if (f != var.samples.end()) {
                    return true;
                } else {
                    return false;
                }
            } else */
            if (type == "Integer") {
                return applyFilter(op, atoi(sample[field].c_str()), atoi(value.c_str()));
            } else if (type == "Float") {
                return applyFilter(op, atof(sample[field].c_str()), atof(value.c_str()));
            } else if (type == "String") {
                return applyFilter(op, sample[field], value);
            } else {
                cerr << "unsupported field type: " << type << endl;
                exit(1);
            }
        }
    }
}

bool VariantCallFile::openVCF(string& filename) {
    open(filename.c_str(), ifstream::in);
    if (!is_open()) {
        cerr << "could not open " << filename << endl;
        return false;
    }
    header = "";
    while (std::getline(*this, line)) {
        if (line.substr(0,2) == "##") {
            // meta-information lines
            // TODO parse into map from info/format key to type
            // ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
            // ##FORMAT=<ID=CB,Number=1,Type=String,Description="Called by S(Sanger), M(UMich), B(BI)">
            size_t found = line.find_first_of("=");
            string entryType = line.substr(2, found - 2);
            // handle reference here, no "<" and ">" given
                //} else if (entryType == "reference") {
            size_t dataStart = line.find_first_of("<");
            size_t dataEnd = line.find_first_of(">");
            if (dataStart != string::npos && dataEnd != string::npos) {
                string entryData = line.substr(dataStart + 1, dataEnd - dataStart - 1);
                // XXX bad; this will break if anyone ever moves the order
                // of the fields around to include a "long form" string
                // including either a = or , in the first or second field
                if (entryType == "INFO" || entryType == "FORMAT") {
                    vector<string> fields = split(entryData, "=,");
                    assert(fields[0] == "ID");
                    string id = fields[1];
                    assert(fields[2] == "Number");
                    int number = atoi(fields[3].c_str());
                    assert(fields[4] == "Type");
                    string type = fields[5];
                    if (entryType == "INFO") {
                        infoCounts[id] = number;
                        infoTypes[id] = type;
                    } else if (entryType == "FORMAT") {
                        formatCounts[id] = number;
                        formatTypes[id] = type;
                    }
                }
            }
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

bool VariantCallFile::getNextVariant(Variant& var) {
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

