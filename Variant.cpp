#include "Variant.h"

namespace vcf {

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
    // if the format changed, we have to rebuild the samples
    if (fields.at(8) != lastFormat) {
        samples.clear();
        lastFormat = fields.at(8);
    }
    vector<string>::iterator sampleName = sampleNames.begin();
    vector<string>::iterator sample = fields.begin() + 9;
    for (; sample != fields.end() && sampleName != sampleNames.end(); ++sample, ++sampleName) {
        string& name = *sampleName;
        if (*sample == ".") {
            samples.erase(name);
            continue;
        }
        vector<string> samplefields = split(*sample, ':');
        vector<string>::iterator i = samplefields.begin();
        for (vector<string>::iterator f = format.begin(); f != format.end(); ++f) {
            samples[name][*f] = *i; ++i;
        }
    }
    if (sampleName != sampleNames.end()) {
        cerr << "error: more sample names in header than sample fields" << endl;
        cerr << "samples: " << join(sampleNames, " ") << endl;
        cerr << "line: " << line << endl;
        exit(1);
    }
    if (sample != fields.end()) {
        cerr << "error: more sample fields than samples listed in header" << endl;
        cerr << "samples: " << join(sampleNames, " ") << endl;
        cerr << "line: " << line << endl;
        exit(1);
    }
    //return true; // we should be catching exceptions...
}

bool Variant::getInfoValueBool(string& key) {
    map<string, string>::iterator s = vcf.infoTypes.find(key);
    if (s == vcf.infoTypes.end()) {
        cerr << "no info field " << key << endl;
        exit(1);
    } else {
        string& type = s->second;
        if (type == "Flag") {
            map<string, bool>::iterator b = infoFlags.find(key);
            if (b == infoFlags.end())
                return false;
            else
                return true;
        } else {
            cerr << "not flag type " << key << endl;
        }
    }
}

string Variant::getInfoValueString(string& key) {
    map<string, string>::iterator s = vcf.infoTypes.find(key);
    if (s == vcf.infoTypes.end()) {
        cerr << "no info field " << key << endl;
        exit(1);
    } else {
        string& type = s->second;
        if (type == "String") {
            map<string, string>::iterator b = info.find(key);
            if (b == info.end())
                return "";
            return b->second;
        } else {
            cerr << "not string type " << key << endl;
            return "";
        }
    }
}

float Variant::getInfoValueFloat(string& key) {
    map<string, string>::iterator s = vcf.infoTypes.find(key);
    if (s == vcf.infoTypes.end()) {
        cerr << "no info field " << key << endl;
        exit(1);
    } else {
        string& type = s->second;
        if (type == "Integer") {
            map<string, string>::iterator b = info.find(key);
            if (b == info.end())
                return false;
            int r;
            if (!convert(b->second, r)) {
                cerr << "could not convert field " << b->second << " to " << type << endl;
                exit(1);
            }
            return r;
        } else if (type == "Float") {
            map<string, string>::iterator b = info.find(key);
            if (b == info.end())
                return false;
            float r;
            if (!convert(b->second, r)) {
                cerr << "could not convert field " << b->second << " to " << type << endl;
                exit(1);
            }
            return r;
        } else {
            cerr << "unsupported type for variant record " << type << endl;
        }
    }
}

bool Variant::getSampleValueBool(string& key, string& sample) {
    map<string, string>::iterator s = vcf.formatTypes.find(key);
    if (s == vcf.infoTypes.end()) {
        cerr << "no info field " << key << endl;
        exit(1);
    } else {
        string& type = s->second;
        map<string, string>& sampleData = samples[sample];
        if (type == "Flag") {
            map<string, string>::iterator b = sampleData.find(key);
            if (b == sampleData.end())
                return false;
            else
                return true;
        } else {
            cerr << "not string type " << key << endl;
        }
    }
}

string Variant::getSampleValueString(string& key, string& sample) {
    map<string, string>::iterator s = vcf.formatTypes.find(key);
    if (s == vcf.infoTypes.end()) {
        cerr << "no info field " << key << endl;
        exit(1);
    } else {
        string& type = s->second;
        map<string, string>& sampleData = samples[sample];
        if (type == "String") {
            map<string, string>::iterator b = sampleData.find(key);
            if (b == sampleData.end())
                return false;
            return b->second;
        } else {
            cerr << "not string type " << key << endl;
        }
    }
}

float Variant::getSampleValueFloat(string& key, string& sample) {
    map<string, string>::iterator s = vcf.formatTypes.find(key);
    if (s == vcf.infoTypes.end()) {
        cerr << "no info field " << key << endl;
        exit(1);
    } else {
        string& type = s->second;
        map<string, string>& sampleData = samples[sample];
        if (type == "Integer") {
            map<string, string>::iterator b = sampleData.find(key);
            if (b == sampleData.end())
                return false;
            int r;
            if (!convert(b->second, r)) {
                cerr << "could not convert field " << b->second << " to " << type << endl;
                exit(1);
            }
            return r;
        } else if (type == "Float") {
            map<string, string>::iterator b = sampleData.find(key);
            if (b == sampleData.end())
                return false;
            float r;
            if (!convert(b->second, r)) {
                cerr << "could not convert field " << b->second << " to " << type << endl;
                exit(1);
            }
            return r;
        } else {
            cerr << "unsupported type for sample " << type << endl;
        }
    }
}

bool Variant::getValueBool(string& key, string& sample) {
    if (sample.length() == 0) { // an empty sample name means
        return getInfoValueBool(key);
    } else {
        return getSampleValueBool(key, sample);
    }
}

float Variant::getValueFloat(string& key, string& sample) {
    if (sample.length() == 0) { // an empty sample name means
        return getInfoValueFloat(key);
    } else {
        return getSampleValueFloat(key, sample);
    }
}

string Variant::getValueString(string& key, string& sample) {
    if (sample.length() == 0) { // an empty sample name means
        return getInfoValueString(key);
    } else {
        return getSampleValueString(key, sample);
    }
}

void Variant::addFilter(string& tag) {
    if (filter == "")
        filter = tag;
    else
        filter += "," + tag;
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



// shunting yard algorithm
void infixToPrefix(queue<RuleToken> tokens, queue<RuleToken>& prefixtokens) {
    stack<RuleToken> ops;
    while (!tokens.empty()) {
        RuleToken& token = tokens.front();
        if (isOperator(token)) {
            //cerr << "found operator " << token.value << endl;
            while (ops.size() > 0 && isOperator(ops.top())
                    && (   (isLeftAssociative(token)  && priority(token) <= priority(ops.top()))
                        || (isRightAssociative(token) && priority(token) <  priority(ops.top())))) {
                prefixtokens.push(ops.top());
                ops.pop();
            }
            ops.push(token);
        } else if (isLeftParenthesis(token)) {
            //cerr << "found paran " << token.value << endl;
            ops.push(token);
        } else if (isRightParenthesis(token)) {
            //cerr << "found paran " << token.value << endl;
            while (ops.size() > 0 && !isLeftParenthesis(ops.top())) {
                prefixtokens.push(ops.top());
                ops.pop();
            }
            if (ops.size() == 0) {
                cerr << "error: mismatched parentheses" << endl;
                exit(1);
            }
            if (isLeftParenthesis(ops.top())) {
                ops.pop();
            }
        } else {
            //cerr << "found operand " << token.value << endl;
            prefixtokens.push(token);
        }
        tokens.pop();
    }
    while (ops.size() > 0) {
        if (isRightParenthesis(ops.top()) || isLeftParenthesis(ops.top())) {
            cerr << "error: mismatched parentheses" << endl;
            exit(1);
        }
        prefixtokens.push(ops.top());
        ops.pop();
    }
}

RuleToken::RuleToken(string tokenstr) {
    isVariable = false;
    if (tokenstr == "!") {
        type = RuleToken::NOT_OPERATOR;
    } else if (tokenstr == "&") {
        type = RuleToken::AND_OPERATOR;
    } else if (tokenstr == "|") {
        type = RuleToken::OR_OPERATOR;
    } else if (tokenstr == "=") {
        type = RuleToken::EQUAL_OPERATOR;
    } else if (tokenstr == ">") {
        type = RuleToken::GREATER_THAN_OPERATOR;
    } else if (tokenstr == "<") {
        type = RuleToken::LESS_THAN_OPERATOR;
    } else if (tokenstr == "(") {
        type = RuleToken::LEFT_PARENTHESIS;
    } else if (tokenstr == ")") {
        type = RuleToken::RIGHT_PARENTHESIS;
    } else { // operand
        type = RuleToken::OPERAND;
        if (convert(tokenstr, number)) {
            type = RuleToken::NUMBER;
        } else {
            // TODO
            isVariable = true;
        }
    }
    value = tokenstr;
}


void tokenizeFilterSpec(string& filterspec, queue<RuleToken>& tokens) {
    string lastToken = "";
    bool inToken = false;
    for (int i = 0; i < filterspec.size(); ++i) {
        char c = filterspec.at(i);
        if (c == ' ') {
            inToken = false;
            if (!inToken && lastToken.size() > 0) {
                tokens.push(RuleToken(lastToken));
                lastToken = "";
            }
        } else if (isOperatorChar(c) || isParanChar(c)) {
            inToken = false;
            if (lastToken.size() > 0) {
                tokens.push(RuleToken(lastToken));
                lastToken = "";
            }
            tokens.push(RuleToken(filterspec.substr(i,1)));
        } else {
            inToken = true;
            lastToken += c;
        }
    }
    // get the last token
    if (inToken) {
        tokens.push(RuleToken(lastToken));
    }
}

// class which evaluates filter expressions
// allow filters to be defined using boolean infix expressions e.g.:
//
// "GQ > 10 & (DP < 3 | DP > 5) & SAMPLE = NA12878"
// or
// "GT = 1/1 | GT = 0/0"
//
// on initialization, tokenizes the input sequence, and converts it from infix to postfix
// on call to 
//


VariantFilter::VariantFilter(string filterspec, VariantFilterType filtertype) {
    type = filtertype;
    spec = filterspec;
    tokenizeFilterSpec(filterspec, tokens);
    infixToPrefix(tokens, rules);
    /*while (!rules.empty()) {
        cerr << " " << rules.front().value << ((isNumeric(rules.front())) ? "f" : "");
        rules.pop();
    }
    */
    //cerr << endl;
    //cerr << join(" ", tokens) << endl;
}

bool VariantFilter::passes(Variant& var, string& sample) {
    // to evaluate a rpn boolean queue with embedded numbers and variables
    // make a result stack, use float to allow comparison of floating point
    // numbers, booleans, and integers
    stack<RuleToken> results;
    queue<RuleToken> rulesCopy = rules; // copy

    while (!rulesCopy.empty()) {
        RuleToken& token = rulesCopy.front();
        rulesCopy.pop();
        // pop operands from the front of the queue and push them onto the stack
        if (isOperand(token)) {
            cerr << "is operand: " << token.value << endl;
            // if the token is variable, i.e. not evaluated in this context, we
            // must evaluate it before pushing it onto the stack
            if (token.isVariable) {
                cerr << "is variable" << endl;
                // look up the variable using the Variant, depending on our filter type
                string type;
                if (sample.empty()) { // means we are record-specific
                    type = var.vcf.infoTypes[token.value];
                } else {
                    type = var.vcf.formatTypes[token.value];
                }
                cerr << "token.value " << token.value << endl;
                cerr << "type: " << type << endl;

                if (type == "Integer" || type == "Float") {
                    token.type = RuleToken::NUMERIC_VARIABLE;
                    token.number = var.getValueFloat(token.value, sample);
                } else if (isString(token)) {
                    token.type = RuleToken::STRING_VARIABLE;
                    token.str = var.getValueString(token.value, sample);
                } else if (type == "Flag") {
                    token.type = RuleToken::BOOLEAN_VARIABLE;
                    token.state = var.getValueBool(token.value, sample);
                }
            }
            results.push(token);
        } 
        // apply operators to the first n elements on the stack and push the result back onto the stack
        else if (isOperator(token)) {
            cerr << "is operator: " << token.value << endl;
            RuleToken a, b, r;
            // is it a not-operator?
            switch (token.type) {
                case ( RuleToken::NOT_OPERATOR ):
                    a = results.top();
                    results.pop();
                    if (!isBoolean(a)) {
                        cerr << "cannot negate a non-boolean" << endl;
                    } else {
                        a.state = !a.state;
                        results.push(a);
                    }
                    break;
                case ( RuleToken::EQUAL_OPERATOR ):
                    a = results.top(); results.pop();
                    b = results.top(); results.pop();
                    if (a.type == b.type) {
                        switch (a.type) {
                            case (RuleToken::STRING_VARIABLE):
                                r.state = (a.str == b.str);
                                break;
                            case (RuleToken::NUMERIC_VARIABLE):
                                r.state = (a.number == b.number);
                                break;
                            case (RuleToken::BOOLEAN_VARIABLE):
                                r.state = (a.state == b.state);
                                break;
                            default:
                                cerr << "should not get here" << endl; exit(1);
                                break;
                        }
                    }
                    results.push(r);
                    break;
                case ( RuleToken::GREATER_THAN_OPERATOR ):
                case ( RuleToken::LESS_THAN_OPERATOR ):
                    a = results.top(); results.pop();
                    b = results.top(); results.pop();
                    if (a.type == b.type && a.type == RuleToken::NUMERIC_VARIABLE) {
                        if (token.type == RuleToken::GREATER_THAN_OPERATOR) {
                            r.state = (a.number > b.number);
                        } else {
                            r.state = (a.number < b.number);
                        }
                    } else {
                        cerr << "cannot compare (> or <) objects of dissimilar types" << endl;
                        exit(1);
                    }
                    results.push(r);
                    break;
                case ( RuleToken::AND_OPERATOR ):
                case ( RuleToken::OR_OPERATOR ):
                    a = results.top(); results.pop();
                    b = results.top(); results.pop();
                    if (a.type == b.type && a.type == RuleToken::BOOLEAN_VARIABLE) {
                        if (token.type == RuleToken::AND_OPERATOR) {
                            r.state = (a.state && b.state);
                        } else {
                            r.state = (a.state || b.state);
                        }
                    } else {
                        cerr << "cannot compare (& or |) objects of dissimilar types" << endl;
                        exit(1);
                    }
                    results.push(r);
                    break;
                default:
                    cerr << "should not get here!" << endl; exit(1);
                    break;
            }
        }
    }
    // at the end you should have only one value on the stack, return it as a boolean
    if (results.size() == 1) {
        if (isBoolean(results.top())) {
           return results.top().state;
        } else {
            cerr << "error, non-boolean value left on stack" << endl;
            //cerr << results.top().value << endl;
            exit(1);
        }
    } else if (results.size() > 1) {
        cerr << "more than one value left on results stack!" << endl;
        while (!results.empty()) {
            cerr << results.top().value << endl;
            results.pop();
        }
        exit(1);
    } else {
        cerr << "results stack empty" << endl;
        exit(1);
    }
}

/*
bool VariantCallFile::openVCF(string& filename) {
    file.open(filename.c_str(), ifstream::in);
    if (!file.is_open()) {
        cerr << "could not open " << filename << endl;
        return false;
    } else {
        return parseHeader();
    }
}

bool VariantCallFile::openVCF(ifstream& stream) {
    file = stream;
    if (!file.is_open()) {
        cerr << "provided file is not open" << endl;
        return false;
    } else {
        return parseHeader();
    }
}
*/

bool VariantCallFile::parseHeader(void) {
    header = "";
    while (std::getline(*file, line)) {
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
                        //cerr << id << " == " << type << endl;
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
    if (std::getline(*file, line)) {
        var.parse(line);
        return true;
    } else {
        return false;
    }
}

} // end namespace vcf
