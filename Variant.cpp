#include "Variant.h"

namespace vcf {

void Variant::parse(string& line, bool parseSamples) {

    // clean up potentially variable data structures
    info.clear();
    infoFlags.clear();
    format.clear();
    alt.clear();
    alleles.clear();

    // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT [SAMPLE1 .. SAMPLEN]
    vector<string> fields = split(line, '\t');

    sequenceName = fields.at(0);
    char* end; // dummy variable for strtoll
    position = strtoll(fields.at(1).c_str(), &end, 10);
    id = fields.at(2);
    ref = fields.at(3);
    alt = split(fields.at(4), ","); // a comma-separated list of alternate alleles

    // make a list of all (ref + alts) alleles, allele[0] = ref, alleles[1:] = alts
    // add the ref allele ([0]), resize for the alt alleles, and then add the alt alleles
    alleles.push_back(ref);
    alleles.resize(alt.size()+1);
    std::copy(alt.begin(), alt.end(), alleles.begin()+1);

    // set up reverse lookup of allele index
    altAlleleIndexes.clear();
    int i = 0;
    for (vector<string>::iterator a = alt.begin();
            a != alt.end(); ++a, ++i) {
        altAlleleIndexes[*a] = i;
    }

    convert(fields.at(5), quality);
    filter = fields.at(6);
    if (fields.size() > 7) {
	vector<string> infofields = split(fields.at(7), ';');
	for (vector<string>::iterator f = infofields.begin(); f != infofields.end(); ++f) {
	    if (*f == ".") {
		continue;
	    }
	    vector<string> kv = split(*f, '=');
	    if (kv.size() == 2) {
		info[kv.at(0)] = split(kv.at(1), ',');
	    } else if (kv.size() == 1) {
		infoFlags[kv.at(0)] = true;
	    }
	}
    }
    // check if we have samples specified
    // and that we are supposed to parse them
    if (parseSamples && fields.size() > 8) {
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
            if (*sample == "." || *sample == "./.") {
                samples.erase(name);
                continue;
            }
            vector<string> samplefields = split(*sample, ':');
            vector<string>::iterator i = samplefields.begin();
            if (samplefields.size() != format.size()) {
                // ignore this case... malformed (or 'null') sample specs are caught above
                // /*
                // cerr << "inconsistent number of fields for sample " << name << endl
                //      << "format is " << join(format, ":") << endl
                //      << "sample is " << *sample << endl;
                // exit(1);
                // *
            }
            else {
                for (vector<string>::iterator f = format.begin(); f != format.end(); ++f) {
                    samples[name][*f] = split(*i, ','); ++i;
                }
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
    }

    //return true; // we should be catching exceptions...
}

void Variant::setVariantCallFile(VariantCallFile& v) {
    sampleNames = v.sampleNames;
    outputSampleNames = v.sampleNames;
    vcf = &v;
}

void Variant::setVariantCallFile(VariantCallFile* v) {
    sampleNames = v->sampleNames;
    outputSampleNames = v->sampleNames;
    vcf = v;
}

ostream& operator<<(ostream& out, VariantFieldType type) {
    switch (type) {
        case FIELD_INTEGER:
            out << "integer";
            break;
        case FIELD_FLOAT:
            out << "float";
            break;
        case FIELD_BOOL:
            out << "bool";
            break;
        case FIELD_STRING:
            out << "string";
            break;
        default:
            out << "unknown";
            break;
    }
    return out;
}

VariantFieldType typeStrToVariantFieldType(string& typeStr) {
    if (typeStr == "Integer") {
        return FIELD_INTEGER;
    } else if (typeStr == "Float") {
        return FIELD_FLOAT;
    } else if (typeStr == "Flag") {
        return FIELD_BOOL;
    } else if (typeStr == "String") {
        return FIELD_STRING;
    } else {
        return FIELD_UNKNOWN;
    }
}

VariantFieldType Variant::infoType(string& key) {
    map<string, VariantFieldType>::iterator s = vcf->infoTypes.find(key);
    if (s == vcf->infoTypes.end()) {
        if (key == "QUAL") { // hack to use QUAL as an "info" field
            return FIELD_INTEGER;
        }
        cerr << "no info field " << key << endl;
        exit(1);
    } else {
        return s->second;
    }
}

VariantFieldType Variant::formatType(string& key) {
    map<string, VariantFieldType>::iterator s = vcf->formatTypes.find(key);
    if (s == vcf->formatTypes.end()) {
        cerr << "no format field " << key << endl;
        exit(1);
    } else {
        return s->second;
    }
}

bool Variant::getInfoValueBool(string& key, int index) {
    map<string, VariantFieldType>::iterator s = vcf->infoTypes.find(key);
    if (s == vcf->infoTypes.end()) {
        cerr << "no info field " << key << endl;
        exit(1);
    } else {
        int count = vcf->infoCounts[key];
        // XXX TODO, fix for Genotype variants...
        if (count != ALLELE_NUMBER) {
            index = 0;
        }
        if (index == INDEX_NONE) {
            if (count != 1) {
                cerr << "no field index supplied and field count != 1" << endl;
                exit(1);
            } else {
                index = 0;
            }
        }
        VariantFieldType type = s->second;
        if (type == FIELD_BOOL) {
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

string Variant::getInfoValueString(string& key, int index) {
    map<string, VariantFieldType>::iterator s = vcf->infoTypes.find(key);
    if (s == vcf->infoTypes.end()) {
        cerr << "no info field " << key << endl;
        exit(1);
    } else {
        int count = vcf->infoCounts[key];
        // XXX TODO, fix for Genotype variants...
        if (count != ALLELE_NUMBER) {
            index = 0;
        }
        if (index == INDEX_NONE) {
            if (count != 1) {
                cerr << "no field index supplied and field count != 1" << endl;
                exit(1);
            } else {
                index = 0;
            }
        }
        VariantFieldType type = s->second;
        if (type == FIELD_STRING) {
            map<string, vector<string> >::iterator b = info.find(key);
            if (b == info.end())
                return "";
            return b->second.at(index);
        } else {
            cerr << "not string type " << key << endl;
            return "";
        }
    }
}

double Variant::getInfoValueFloat(string& key, int index) {
    map<string, VariantFieldType>::iterator s = vcf->infoTypes.find(key);
    if (s == vcf->infoTypes.end()) {
        if (key == "QUAL") {
            return quality;
        }
        cerr << "no info field " << key << endl;
        exit(1);
    } else {
        int count = vcf->infoCounts[key];
        // XXX TODO, fix for Genotype variants...
        if (count != ALLELE_NUMBER) {
            index = 0;
        }
        if (index == INDEX_NONE) {
            if (count != 1) {
                cerr << "no field index supplied and field count != 1" << endl;
                exit(1);
            } else {
                index = 0;
            }
        }
        VariantFieldType type = s->second;
        if (type == FIELD_FLOAT || type == FIELD_INTEGER) {
            map<string, vector<string> >::iterator b = info.find(key);
            if (b == info.end())
                return false;
            double r;
            if (!convert(b->second.at(index), r)) {
                cerr << "could not convert field " << key << "=" << b->second.at(index) << " to " << type << endl;
                exit(1);
            }
            return r;
        } else {
            cerr << "unsupported type for variant record " << type << endl;
        }
    }
}

int Variant::getNumSamples(void) {
    return sampleNames.size();
}

int Variant::getNumValidGenotypes(void) {
    int valid_genotypes = 0;
    map<string, map<string, vector<string> > >::const_iterator s     = samples.begin();
    map<string, map<string, vector<string> > >::const_iterator sEnd  = samples.end();
    for (; s != sEnd; ++s) {
        map<string, vector<string> > sample_info = s->second;
        if (sample_info["GT"].front() != "./.") {
            valid_genotypes++;
        }
    }
    return valid_genotypes;
}

bool Variant::getSampleValueBool(string& key, string& sample, int index) {
    map<string, VariantFieldType>::iterator s = vcf->formatTypes.find(key);
    if (s == vcf->infoTypes.end()) {
        cerr << "no info field " << key << endl;
        exit(1);
    } else {
        int count = vcf->formatCounts[key];
        // XXX TODO, fix for Genotype variants...
        if (count != ALLELE_NUMBER) {
            index = 0;
        }
        if (index == INDEX_NONE) {
            if (count != 1) {
                cerr << "no field index supplied and field count != 1" << endl;
                exit(1);
            } else {
                index = 0;
            }
        }
        VariantFieldType type = s->second;
        map<string, vector<string> >& sampleData = samples[sample];
        if (type == FIELD_BOOL) {
            map<string, vector<string> >::iterator b = sampleData.find(key);
            if (b == sampleData.end())
                return false;
            else
                return true;
        } else {
            cerr << "not bool type " << key << endl;
        }
    }
}

string Variant::getSampleValueString(string& key, string& sample, int index) {
    map<string, VariantFieldType>::iterator s = vcf->formatTypes.find(key);
    if (s == vcf->infoTypes.end()) {
        cerr << "no info field " << key << endl;
        exit(1);
    } else {
        int count = vcf->formatCounts[key];
        // XXX TODO, fix for Genotype variants...
        if (count != ALLELE_NUMBER) {
            index = 0;
        }
        if (index == INDEX_NONE) {
            if (count != 1) {
                cerr << "no field index supplied and field count != 1" << endl;
                exit(1);
            } else {
                index = 0;
            }
        }
        VariantFieldType type = s->second;
        map<string, vector<string> >& sampleData = samples[sample];
        if (type == FIELD_STRING) {
            map<string, vector<string> >::iterator b = sampleData.find(key);
            if (b == sampleData.end()) {
                return "";
            } else {
                return b->second.at(index);
            }
        } else {
            cerr << "not string type " << key << endl;
        }
    }
}

double Variant::getSampleValueFloat(string& key, string& sample, int index) {
    map<string, VariantFieldType>::iterator s = vcf->formatTypes.find(key);
    if (s == vcf->infoTypes.end()) {
        cerr << "no info field " << key << endl;
        exit(1);
    } else {
        // XXX TODO wrap this with a function call
        int count = vcf->formatCounts[key];
        // XXX TODO, fix for Genotype variants...
        if (count != ALLELE_NUMBER) {
            index = 0;
        }
        if (index == INDEX_NONE) {
            if (count != 1) {
                cerr << "no field index supplied and field count != 1" << endl;
                exit(1);
            } else {
                index = 0;
            }
        }
        VariantFieldType type = s->second;
        map<string, vector<string> >& sampleData = samples[sample];
        if (type == FIELD_FLOAT || type == FIELD_INTEGER) {
            map<string, vector<string> >::iterator b = sampleData.find(key);
            if (b == sampleData.end())
                return false;
            double r;
            if (!convert(b->second.at(index), r)) {
                cerr << "could not convert field " << key << "=" << b->second.at(index) << " to " << type << endl;
                exit(1);
            }
            return r;
        } else {
            cerr << "unsupported type for sample " << type << endl;
        }
    }
}

bool Variant::getValueBool(string& key, string& sample, int index) {
    if (sample.empty()) { // an empty sample name means
        return getInfoValueBool(key, index);
    } else {
        return getSampleValueBool(key, sample, index);
    }
}

double Variant::getValueFloat(string& key, string& sample, int index) {
    if (sample.empty()) { // an empty sample name means
        return getInfoValueFloat(key, index);
    } else {
        return getSampleValueFloat(key, sample, index);
    }
}

string Variant::getValueString(string& key, string& sample, int index) {
    if (sample.empty()) { // an empty sample name means
        return getInfoValueString(key, index);
    } else {
        return getSampleValueString(key, sample, index);
    }
}

int Variant::getAltAlleleIndex(string& allele) {
    map<string, int>::iterator f = altAlleleIndexes.find(allele);
    if (f == altAlleleIndexes.end()) {
        cerr << "no such allele \'" << allele << "\' in record " << sequenceName << ":" << position << endl;
        exit(1);
    } else {
        return f->second;
    }
}

void Variant::addFilter(string& tag) {
    if (filter == "" || filter == ".")
        filter = tag;
    else
        filter += "," + tag;
}

void Variant::addFormatField(string& key) {
    bool hasTag = false;
    for (vector<string>::iterator t = format.begin(); t != format.end(); ++t) {
        if (*t == key) {
            hasTag = true;
            break;
        }
    }
    if (!hasTag) {
        format.push_back(key);
    }
}

void Variant::printAlt(ostream& out) {
    for (vector<string>::iterator i = alt.begin(); i != alt.end(); ++i) {
        out << *i;
        // add a comma for all but the last alternate allele
        if (i != (alt.end() - 1)) out << ",";
    }
}

void Variant::printAlleles(ostream& out) {
    for (vector<string>::iterator i = alleles.begin(); i != alleles.end(); ++i) {
        out << *i;
        // add a comma for all but the last alternate allele
        if (i != (alleles.end() - 1)) out << ",";
    }
}

ostream& operator<<(ostream& out, Variant& var) {
    out << var.sequenceName << "\t"
        << var.position << "\t"
        << var.id << "\t"
        << var.ref << "\t";
    // report the list of alternate alleles.
    var.printAlt(out);
    out << "\t"
        << var.quality << "\t"
        << var.filter << "\t";
    for (map<string, vector<string> >::iterator i = var.info.begin(); i != var.info.end(); ++i) {
        out << ((i == var.info.begin()) ? "" : ";") << i->first << "=" << join(i->second, ",");
    }
    for (map<string, bool>::iterator i = var.infoFlags.begin(); i != var.infoFlags.end(); ++i) {
        if (i == var.infoFlags.end()) {
            out << "";
        } else if (i == var.infoFlags.begin() && var.info.empty()) {
            out << "";
        } else {
            out << ";";
        }
        out << i->first;
    }
    if (!var.format.empty()) {
        out << "\t";
        for (vector<string>::iterator f = var.format.begin(); f != var.format.end(); ++f) {
            out << ((f == var.format.begin()) ? "" : ":") << *f;
        }
        for (vector<string>::iterator s = var.outputSampleNames.begin(); s != var.outputSampleNames.end(); ++s) {
            out << "\t";
            map<string, map<string, vector<string> > >::iterator sampleItr = var.samples.find(*s);
            if (sampleItr == var.samples.end()) {
                out << ".";
            } else {
                map<string, vector<string> >& sample = sampleItr->second;
                if (sample.size() == 0) {
                    out << ".";
                } else {
                    for (vector<string>::iterator f = var.format.begin(); f != var.format.end(); ++f) {
                        map<string, vector<string> >::iterator g = sample.find(*f);
                        out << ((f == var.format.begin()) ? "" : ":");
                        if (g != sample.end()) {
                            out << join(g->second, ",");
                        } else {
                            out << ".";
                        }
                    }
                }
            }
        }
    }
    return out;
}

void Variant::setOutputSampleNames(vector<string>& samplesToOutput) {
    outputSampleNames = samplesToOutput;
}

float Variant::getAAF(void) {
    
    if (alleles.size() > 2)
        return -1.0;
        
    int num_valid_genotypes = 0;
    int tot_alt_alleles = 0;
    map<string, map<string, vector<string> > >::const_iterator s     = samples.begin();
    map<string, map<string, vector<string> > >::const_iterator sEnd  = samples.end();
    for (; s != sEnd; ++s) {
        map<string, vector<string> > sample_info = s->second;
        int gt_val = genotypeValue(sample_info["GT"].front());
        if (gt_val >= 0) {
            num_valid_genotypes++;
            tot_alt_alleles += gt_val;
        }
    }
    return ((float) tot_alt_alleles / (float) (2.0 * num_valid_genotypes));
}

float Variant::getNucleotideDiversity(void) {
    
    int num_valid_genotypes = 0;
    int tot_alt_alleles = 0;
    map<string, map<string, vector<string> > >::const_iterator s     = samples.begin();
    map<string, map<string, vector<string> > >::const_iterator sEnd  = samples.end();
    for (; s != sEnd; ++s) {
        map<string, vector<string> > sample_info = s->second;
        int gt_val = genotypeValue(sample_info["GT"].front());
        if (gt_val >= 0) {
            num_valid_genotypes++;
            tot_alt_alleles += gt_val;
        }
    }
    int num_chroms = 2 * num_valid_genotypes;
    float p = ((float) tot_alt_alleles / (float) (2.0 * num_chroms));
    float q = 1.0 - p;
    return (float) (num_chroms/(num_chroms-1.0)) * (2.0 * p * q);
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

RuleToken::RuleToken(string tokenstr, map<string, VariantFieldType>& variables) {
    isVariable = false;
    if (tokenstr == "!") {
        type = RuleToken::NOT_OPERATOR;
    } else if (tokenstr == "&") {
        type = RuleToken::AND_OPERATOR;
    } else if (tokenstr == "|") {
        type = RuleToken::OR_OPERATOR;
    } else if (tokenstr == "+") {
        type = RuleToken::ADD_OPERATOR;
    } else if (tokenstr == "-") {
        type = RuleToken::SUBTRACT_OPERATOR;
    } else if (tokenstr == "*") {
        type = RuleToken::MULTIPLY_OPERATOR;
    } else if (tokenstr == "/") {
        type = RuleToken::DIVIDE_OPERATOR;
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
        if (variables.find(tokenstr) == variables.end()) {
            if (convert(tokenstr, number)) {
                type = RuleToken::NUMBER;
            } else if (tokenstr == "QUAL") {
                isVariable = true;
            } else {
                type = RuleToken::STRING_VARIABLE;
            }
        } else {
            isVariable = true;
        }
    }
    value = tokenstr;
}


void tokenizeFilterSpec(string& filterspec, queue<RuleToken>& tokens, map<string, VariantFieldType>& variables) {
    string lastToken = "";
    bool inToken = false;
    for (int i = 0; i < filterspec.size(); ++i) {
        char c = filterspec.at(i);
        if (c == ' ' || c == '\n') {
            inToken = false;
            if (!inToken && lastToken.size() > 0) {
                tokens.push(RuleToken(lastToken, variables));
                lastToken = "";
            }
        } else if (!inToken && (isOperatorChar(c) || isParanChar(c))) {
            inToken = false;
            if (lastToken.size() > 0) {
                tokens.push(RuleToken(lastToken, variables));
                lastToken = "";
            }
            tokens.push(RuleToken(filterspec.substr(i,1), variables));
        } else {
            inToken = true;
            lastToken += c;
        }
    }
    // get the last token
    if (inToken) {
        tokens.push(RuleToken(lastToken, variables));
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


VariantFilter::VariantFilter(string filterspec, VariantFilterType filtertype, map<string, VariantFieldType>& variables) {
    type = filtertype;
    spec = filterspec;
    tokenizeFilterSpec(filterspec, tokens, variables);
    infixToPrefix(tokens, rules);
    /*while (!rules.empty()) {
        cerr << " " << rules.front().value << ((isNumeric(rules.front())) ? "f" : "");
        rules.pop();
    }
    */
    //cerr << endl;
    //cerr << join(" ", tokens) << endl;
}

// all alts pass
bool VariantFilter::passes(Variant& var, string& sample) {
    for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
        string& allele = *a;
        if (!passes(var, sample, allele)) {
            return false;
        }
    }
    return true;
}

bool VariantFilter::passes(Variant& var, string& sample, string& allele) {
    // to evaluate a rpn boolean queue with embedded numbers and variables
    // make a result stack, use float to allow comparison of floating point
    // numbers, booleans, and integers
    stack<RuleToken> results;
    queue<RuleToken> rulesCopy = rules; // copy

    int index;
    if (allele.empty()) {
        index = 0; // apply to the whole record
    } else {
        // apply to a specific allele
        index = var.getAltAlleleIndex(allele);
    }

    while (!rulesCopy.empty()) {
        RuleToken& token = rulesCopy.front();
        rulesCopy.pop();
        // pop operands from the front of the queue and push them onto the stack
        if (isOperand(token)) {
            //cout << "is operand: " << token.value << endl;
            // if the token is variable, i.e. not evaluated in this context, we
            // must evaluate it before pushing it onto the stack
            if (token.isVariable) {
                //cout << "is variable" << endl;
                // look up the variable using the Variant, depending on our filter type
                //cout << "token.value " << token.value << endl;
                VariantFieldType type;
                if (sample.empty()) { // means we are record-specific
                    type = var.infoType(token.value);
                } else {
                    type = var.formatType(token.value);
                    //cout << "type = " << type << endl;
                }
                //cout << "type: " << type << endl;

                if (type == FIELD_INTEGER || type == FIELD_FLOAT) {
                    token.type = RuleToken::NUMERIC_VARIABLE;
                    token.number = var.getValueFloat(token.value, sample, index);
                    //cerr << "number: " << token.number << endl;
                } else if (type == FIELD_BOOL) {
                    token.type = RuleToken::BOOLEAN_VARIABLE;
                    token.state = var.getValueBool(token.value, sample, index);
                    //cerr << "state: " << token.state << endl;
                } else if (type == FIELD_STRING) {
                    //cout << "token.value = " << token.value << endl;
                    token.type = RuleToken::STRING_VARIABLE;
                    token.str = var.getValueString(token.value, sample, index);
                } else if (isString(token)) {
                    token.type = RuleToken::STRING_VARIABLE;
                    token.str = var.getValueString(token.value, sample, index);
                    //cerr << "string: " << token.str << endl;
                }
            } else {
                double f;
                string s;
                //cerr << "parsing operand" << endl;
                if (convert(token.value, f)) {
                    token.type = RuleToken::NUMERIC_VARIABLE;
                    token.number = f;
                    //cerr << "number: " << token.number << endl;
                } else if (convert(token.value, s)) {
                    token.type = RuleToken::STRING_VARIABLE;
                    token.str = s;
                    //cerr << "string: " << token.str << endl;
                } else {
                    cerr << "could not parse non-variable operand " << token.value << endl;
                    exit(1);
                }
            }
            results.push(token);
        } 
        // apply operators to the first n elements on the stack and push the result back onto the stack
        else if (isOperator(token)) {
            //cerr << "is operator: " << token.value << endl;
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
                    a = results.top(); results.pop();
                    b = results.top(); results.pop();
                    if (a.type == b.type && a.type == RuleToken::NUMERIC_VARIABLE) {
                        r.state = (b.number > a.number);
                    } else {
                        cerr << "cannot compare (>) objects of dissimilar types" << endl;
                        cerr << a.type << " " << b.type << endl;
                        exit(1);
                    }
                    results.push(r);
                    break;

                case ( RuleToken::LESS_THAN_OPERATOR ):
                    a = results.top(); results.pop();
                    b = results.top(); results.pop();
                    if (a.type == b.type && a.type == RuleToken::NUMERIC_VARIABLE) {
                        r.state = (b.number < a.number);
                    } else {
                        cerr << "cannot compare (<) objects of dissimilar types" << endl;
                        cerr << a.type << " " << b.type << endl;
                        exit(1);
                    }
                    results.push(r);
                    break;

                case ( RuleToken::ADD_OPERATOR ):
                    a = results.top(); results.pop();
                    b = results.top(); results.pop();
                    if (a.type == b.type && a.type == RuleToken::NUMERIC_VARIABLE) {
                        r.number = (b.number + a.number);
                        r.type = RuleToken::NUMERIC_VARIABLE;
                    } else {
                        cerr << "cannot add objects of dissimilar types" << endl;
                        cerr << a.type << " " << b.type << endl;
                        exit(1);
                    }
                    results.push(r);
                    break;

                case ( RuleToken::SUBTRACT_OPERATOR ):
                    a = results.top(); results.pop();
                    b = results.top(); results.pop();
                    if (a.type == b.type && a.type == RuleToken::NUMERIC_VARIABLE) {
                        r.number = (b.number - a.number);
                        r.type = RuleToken::NUMERIC_VARIABLE;
                    } else {
                        cerr << "cannot subtract objects of dissimilar types" << endl;
                        cerr << a.type << " " << b.type << endl;
                        exit(1);
                    }
                    results.push(r);
                    break;

                case ( RuleToken::MULTIPLY_OPERATOR ):
                    a = results.top(); results.pop();
                    b = results.top(); results.pop();
                    if (a.type == b.type && a.type == RuleToken::NUMERIC_VARIABLE) {
                        r.number = (b.number * a.number);
                        r.type = RuleToken::NUMERIC_VARIABLE;
                    } else {
                        cerr << "cannot multiply objects of dissimilar types" << endl;
                        cerr << a.type << " " << b.type << endl;
                        exit(1);
                    }
                    results.push(r);
                    break;

                case ( RuleToken::DIVIDE_OPERATOR):
                    a = results.top(); results.pop();
                    b = results.top(); results.pop();
                    if (a.type == b.type && a.type == RuleToken::NUMERIC_VARIABLE) {
                        r.number = (b.number / a.number);
                        r.type = RuleToken::NUMERIC_VARIABLE;
                    } else {
                        cerr << "cannot divide objects of dissimilar types" << endl;
                        cerr << a.type << " " << b.type << endl;
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

void VariantFilter::removeFilteredGenotypes(Variant& var) {

    for (vector<string>::iterator s = var.sampleNames.begin(); s != var.sampleNames.end(); ++s) {
        string& name = *s;
        if (!passes(var, name)) {
            var.samples.erase(name);
        }
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

void VariantCallFile::updateSamples(vector<string>& newSamples) {
    sampleNames = newSamples;
    // regenerate the last line of the header
    vector<string> headerLines = split(header, '\n');
    vector<string> colnames = split(headerLines.at(headerLines.size() - 1), '\t'); // get the last, update the samples
    vector<string> newcolnames;
    newcolnames.resize(9 + sampleNames.size());
    copy(colnames.begin(), colnames.begin() + 9, newcolnames.begin());
    copy(sampleNames.begin(), sampleNames.end(), newcolnames.begin() + 9);
    headerLines.at(headerLines.size() - 1) = join(newcolnames, "\t");
    header = join(headerLines, "\n");
}

// TODO cleanup, store header lines instead of bulk header
void VariantCallFile::addHeaderLine(string line) {
    vector<string> headerLines = split(header, '\n');
    headerLines.insert(headerLines.end() - 1, line);
    header = join(unique(headerLines), "\n");
}

// helper to addHeaderLine
vector<string>& unique(vector<string>& strings) {
    set<string> uniq;
    vector<string> res;
    for (vector<string>::const_iterator s = strings.begin(); s != strings.end(); ++s) {
	if (uniq.find(*s) == uniq.end()) {
	    res.push_back(*s);
	    uniq.insert(*s);
	}
    }
    strings = res;
    return strings;
}

vector<string> VariantCallFile::infoIds(void) {
    vector<string> tags;
    vector<string> headerLines = split(header, '\n');
    for (vector<string>::iterator s = headerLines.begin(); s != headerLines.end(); ++s) {
        string& line = *s;
        if (line.find("##INFO") == 0) {
	    size_t pos = line.find("ID=");
	    if (pos != string::npos) {
		pos += 3;
		size_t tagend = line.find(",", pos);
		if (tagend != string::npos) {
		    tags.push_back(line.substr(pos, tagend - pos));
		}
	    }
	}
    }
    return tags;
}

vector<string> VariantCallFile::formatIds(void) {
    vector<string> tags;
    vector<string> headerLines = split(header, '\n');
    for (vector<string>::iterator s = headerLines.begin(); s != headerLines.end(); ++s) {
        string& line = *s;
        if (line.find("##FORMAT") == 0) {
	    size_t pos = line.find("ID=");
	    if (pos != string::npos) {
		pos += 3;
		size_t tagend = line.find(",", pos);
		if (tagend != string::npos) {
		    tags.push_back(line.substr(pos, tagend - pos));
		}
	    }
	}
    }
    return tags;
}

void VariantCallFile::removeInfoHeaderLine(string tag) {
    vector<string> headerLines = split(header, '\n');
    vector<string> newHeader;
    string id = "ID=" + tag;
    for (vector<string>::iterator s = headerLines.begin(); s != headerLines.end(); ++s) {
        string& line = *s;
        if (line.find("##INFO") == 0) {
            if (line.find(id) == string::npos) {
                newHeader.push_back(line);
            }
        } else {
            newHeader.push_back(line);
        }
    }
    header = join(newHeader, "\n");
}

void VariantCallFile::removeGenoHeaderLine(string tag) {
    vector<string> headerLines = split(header, '\n');
    vector<string> newHeader;
    string id = "ID=" + tag;
    for (vector<string>::iterator s = headerLines.begin(); s != headerLines.end(); ++s) {
        string& line = *s;
        if (line.find("##FORMAT") == 0) {
            if (line.find(id) == string::npos) {
                newHeader.push_back(line);
            }
        } else {
            newHeader.push_back(line);
        }
    }
    header = join(newHeader, "\n");
}

bool VariantCallFile::parseHeader(void) {

    string headerStr = "";

    if (usingTabix) {
        tabixFile->getHeader(headerStr);
        if (headerStr.empty()) {
            cerr << "error: no VCF header" << endl;
            exit(1);
        }
        tabixFile->getNextLine(line);
        firstRecord = true;
    } else {
        while (std::getline(*file, line)) {
            if (line.substr(0,1) == "#") {
                headerStr += line + '\n';
            } else {
                // done with header
                if (headerStr.empty()) {
                    cerr << "error: no VCF header" << endl;
                    exit(1);
                }
                firstRecord = true;
                break;
            }
        }
    }

    return parseHeader(headerStr);

}

bool VariantCallFile::parseHeader(string& h) {

    header = h; // stores the header in the object instance

    vector<string> headerLines = split(header, "\n");
    for (vector<string>::iterator h = headerLines.begin(); h != headerLines.end(); ++h) {
        string headerLine = *h;
        if (headerLine.substr(0,2) == "##") {
            // meta-information headerLines
            // TODO parse into map from info/format key to type
            // ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
            // ##FORMAT=<ID=CB,Number=1,Type=String,Description="Called by S(Sanger), M(UMich), B(BI)">
            size_t found = headerLine.find_first_of("=");
            string entryType = headerLine.substr(2, found - 2);
            // handle reference here, no "<" and ">" given
                //} else if (entryType == "reference") {
            size_t dataStart = headerLine.find_first_of("<");
            size_t dataEnd = headerLine.find_first_of(">");
            if (dataStart != string::npos && dataEnd != string::npos) {
                string entryData = headerLine.substr(dataStart + 1, dataEnd - dataStart - 1);
                // XXX bad; this will break if anyone ever moves the order
                // of the fields around to include a "long form" string
                // including either a = or , in the first or second field
                if (entryType == "INFO" || entryType == "FORMAT") {
                    vector<string> fields = split(entryData, "=,");
                    if (fields[0] != "ID") {
                        cerr << "header parse error at:" << endl
                             << "fields[0] != \"ID\"" << endl
                             << headerLine << endl;
                        exit(1);
                    }
                    string id = fields[1];
                    if (fields[2] != "Number") {
                        cerr << "header parse error at:" << endl
                             << "fields[2] != \"Number\"" << endl
                             << headerLine << endl;
                        exit(1);
                    }
                    int number;
                    string numberstr = fields[3].c_str();
                    // XXX TODO VCF has variable numbers of fields...
                    if (numberstr == "A") {
                        number = ALLELE_NUMBER;
                    } else if (numberstr == "G") {
                        number = GENOTYPE_NUMBER;
		    } else if (numberstr == ".") {
			number = 1;
                    } else {
                        convert(numberstr, number);
                    }
                    if (fields[4] != "Type") {
                        cerr << "header parse error at:" << endl
                             << "fields[4] != \"Type\"" << endl
                             << headerLine << endl;
                        exit(1);
                    }
                    VariantFieldType type = typeStrToVariantFieldType(fields[5]);
                    if (entryType == "INFO") {
                        infoCounts[id] = number;
                        infoTypes[id] = type;
                        //cerr << id << " == " << type << endl;
                    } else if (entryType == "FORMAT") {
                        //cout << "found format field " << id << " with type " << type << endl;
                        formatCounts[id] = number;
                        formatTypes[id] = type;
                    }
                }
            }
        } else if (headerLine.substr(0,1) == "#") {
            // field name headerLine
            vector<string> fields = split(headerLine, '\t');
            if (fields.size() > 8) {
                sampleNames.resize(fields.size() - 9);
                copy(fields.begin() + 9, fields.end(), sampleNames.begin());
            }
        }
    }

    return true;
}

bool VariantCallFile::getNextVariant(Variant& var) {
    if (firstRecord) {
        var.parse(line, parseSamples);
        firstRecord = false;
        _done = false;
        return true;
    }
    if (usingTabix) {
        if (tabixFile->getNextLine(line)) {
            var.parse(line, parseSamples);
            _done = false;
            return true;
        } else {
            _done = true;
            return false;
        }
    } else {
        if (std::getline(*file, line)) {
            var.parse(line, parseSamples);
            _done = false;
            return true;
        } else {
            _done = true;
            return false;
        }
    }
}

bool VariantCallFile::setRegion(string seq, long int start, long int end) {
    stringstream regionstr;
    if (end) {
        regionstr << seq << ":" << start << "-" << end;
    } else {
        regionstr << seq << ":" << start;
    }
    setRegion(regionstr.str());
}

bool VariantCallFile::setRegion(string region) {
    if (!usingTabix) {
        cerr << "cannot setRegion on a non-tabix indexed file" << endl;
        exit(1);
    }
    size_t dots = region.find("..");
    // convert between bamtools/freebayes style region string and tabix/samtools style
    if (dots != string::npos) {
        region.replace(dots, 2, "-");
    }
    if (tabixFile->setRegion(region)) {
        if (tabixFile->getNextLine(line)) {
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}


// genotype manipulation

/*
map<string, int> decomposeGenotype(string& genotype) {
    string splitter = "/";
    if (genotype.find("|") != string::npos) {
        splitter = "|";
    }
    vector<string> haps = split(genotype, splitter);
    map<string, int> decomposed;
    for (vector<string>::iterator h = haps.begin(); h != haps.end(); ++h) {
        ++decomposed[*h];
    }
    return decomposed;
}
*/

map<int, int> decomposeGenotype(string& genotype) {
    string splitter = "/";
    if (genotype.find("|") != string::npos) {
        splitter = "|";
    }
    vector<string> haps = split(genotype, splitter);
    map<int, int> decomposed;
    for (vector<string>::iterator h = haps.begin(); h != haps.end(); ++h) {
        int alt;
        if (*h == ".") {
            ++decomposed[NULL_ALLELE];
        } else {
            convert(*h, alt);
            ++decomposed[alt];
        }
    }
    return decomposed;
}

string genotypeToString(map<int, int>& genotype) {
    vector<int> s;
    for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
        int a = g->first;
        int c = g->second;
        for (int i = 0; i < c; ++i) s.push_back(a);
    }
    sort(s.begin(), s.end());
    vector<string> r;
    for (vector<int>::iterator i = s.begin(); i != s.end(); ++i) {
        if (*i == NULL_ALLELE) r.push_back(".");
        else r.push_back(convert(*i));
    }
    return join(r, "/"); // TODO adjust for phased/unphased
}

bool isHet(map<int, int>& genotype) {
    return genotype.size() > 1;
}

bool isHom(map<int, int>& genotype) {
    return genotype.size() == 1;
}

bool hasNonRef(map<int, int>& genotype) {
    for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
        if (g->first != 0) {
            return true;
        }
    }
    return false;
}

bool isHomRef(map<int, int>& genotype) {
    return isHom(genotype) && !hasNonRef(genotype);
}

bool isHomNonRef(map<int, int>& genotype) {
    return isHom(genotype) && hasNonRef(genotype);
}

bool isNull(map<int, int>& genotype) {
    return genotype.find(NULL_ALLELE) != genotype.end();
}

int ploidy(map<int, int>& genotype) {
    int i = 0;
    for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
        i += g->second;
    }
    return i;
}

int genotypeValue(string & genotype) {
    map<int, int> gt_map = decomposeGenotype(genotype);
    if      (isHomRef(gt_map))    return 0;
    else if (isHet(gt_map))       return 1;
    else if (isHomNonRef(gt_map)) return 2;
    else if (isNull(gt_map))      return -1;
    else cerr << "error, undefined genotype class." << endl;
}

map<string, vector<VariantAllele> > Variant::parsedAlternates(bool includePreviousBaseForIndels) {

    map<string, vector<VariantAllele> > variantAlleles;

    // add the reference allele
    variantAlleles[ref].push_back(VariantAllele(ref, ref, position));

    // single SNP case, no ambiguity possible, no need to spend a lot of
    // compute aligning ref and alt fields
    if (alt.size() == 1 && alt.front().size() == 1) {
        variantAlleles[alt.front()].push_back(VariantAllele(ref, alt.front(), position));
        return variantAlleles;
    }

    // padding is used to ensure a stable alignment of the alternates to the reference
    // without having to go back and look at the full reference sequence
    int paddingLen = max(10, (int) (2 * ref.size()));  // dynamically determine optimum padding length
    for (vector<string>::iterator a = alt.begin(); a != alt.end(); ++a) {
        string& alternate = *a;
	paddingLen = max(paddingLen, (int) (2 * alternate.size()));
    }
    char padChar = 'Z';
    char anchorChar = 'Q';
    string padding(paddingLen, padChar);
    string reference = padding + ref + padding;

    // this 'anchored' string is done for stability
    // the assumption is that there should be a positional match in the first base
    // this is true for VCF 4.1, and standard best practices
    // using the anchor char ensures this without other kinds of realignment
    string reference_M = reference;
    reference_M[paddingLen] = anchorChar;
    //cout << reference_M << endl;

    // passed to sw.Align
    unsigned int referencePos;

    // constants for SmithWaterman algorithm
    float matchScore = 10.0f;
    float mismatchScore = -9.0f;
    float gapOpenPenalty = 15.0f;
    float gapExtendPenalty = 6.66f;

    string cigar;

    for (vector<string>::iterator a = alt.begin(); a != alt.end(); ++a) {

        string& alternate = *a;
        vector<VariantAllele>& variants = variantAlleles[alternate];
        string alternateQuery = padding + alternate + padding;

        // again, the _M string is for stability of alignment aganist VCF alts
        string alternateQuery_M = alternateQuery;
        alternateQuery_M[paddingLen] = anchorChar;
        //const unsigned int alternateLen = alternate.size();

        CSmithWatermanGotoh sw(matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
        sw.Align(referencePos, cigar, reference_M, alternateQuery_M);

        // left-realign the alignment...

        int altpos = 0;
        int refpos = 0;
        int len;
        string slen;
        vector<pair<int, char> > cigarData;

        for (string::iterator c = cigar.begin(); c != cigar.end(); ++c) {
            switch (*c) {
                case 'I':
                    len = atoi(slen.c_str());
                    slen.clear();
                    cigarData.push_back(make_pair(len, *c));
                    if (includePreviousBaseForIndels) {
                        variants.push_back(VariantAllele(reference.substr(refpos - 1, 1), alternateQuery.substr(altpos - 1, len + 1), refpos - paddingLen + position - 1));
                    } else {
                        variants.push_back(VariantAllele("", alternateQuery.substr(altpos, len), refpos - paddingLen + position));
                    }
                    altpos += len;
                    break;
                case 'D':
                    len = atoi(slen.c_str());
                    slen.clear();
                    cigarData.push_back(make_pair(len, *c));
                    if (includePreviousBaseForIndels) {
                        variants.push_back(VariantAllele(reference.substr(refpos - 1, len + 1), alternateQuery.substr(altpos - 1, 1), altpos - paddingLen + position - 1));
                    } else {
                        variants.push_back(VariantAllele(reference.substr(refpos, len), "", refpos - paddingLen + position));
                    }
                    refpos += len;
                    break;
                case 'M':
                    {
                        len = atoi(slen.c_str());
                        slen.clear();
                        cigarData.push_back(make_pair(len, *c));
                        string refmatch = reference.substr(refpos, len);
                        string altmatch = alternateQuery.substr(altpos, len);
                        bool inmismatch = false;
			int mismatchRefPosStart = refpos;
                        int mismatchStart = 0;
                        for (int i = 0; i < refmatch.size(); ++i) {
                            if (refmatch.at(i) == altmatch.at(i)) {
                                if (inmismatch) {
                                    variants.push_back(VariantAllele(
							   refmatch.substr(mismatchStart, i - mismatchStart),
							   altmatch.substr(mismatchStart, i - mismatchStart),
							   mismatchRefPosStart - paddingLen + position));
                                }
                                inmismatch = false;
                            } else {
                                if (!inmismatch) {
				    mismatchRefPosStart = refpos;
				    mismatchStart = i;
                                    inmismatch = true;
                                }
                            }
                            ++refpos;
                            ++altpos;
                        }
			if (inmismatch) {
			    variants.push_back(VariantAllele(
						   refmatch.substr(mismatchStart, refmatch.size() - mismatchStart),
						   altmatch.substr(mismatchStart, refmatch.size() - mismatchStart),
						   mismatchRefPosStart - paddingLen + position));
			}

                    }
                    break;
                case 'S':
                    len = atoi(slen.c_str());
                    slen.clear();
                    cigarData.push_back(make_pair(len, *c));
                    refpos += len;
                    altpos += len;
                    break;
                default:
                    len = 0;
                    slen += *c;
                    break;
            }

	    // if the last two variants have the same alt position,
	    // take the one which covers more ref or alt sequence
	    // this deals with things like ACG/TTCG, which decomposes to A/T, A/TT
	    if (variants.size() > 1) {
		VariantAllele& varA = variants.at(variants.size() - 2);
		VariantAllele& varB = variants.back();
		if (varA.position == varB.position) {
		    if (varA.ref.size() == varB.ref.size()) {
			if (varA.alt.size() >= varB.alt.size()) {
			    variants.pop_back();
			} else {
			    VariantAllele varB_copy = variants.back();
			    variants.pop_back(); variants.pop_back();
			    variants.push_back(varB_copy);
			}
		    } else if (varA.ref.size() > varB.ref.size()) {
			variants.pop_back();
		    } else {
			VariantAllele varB_copy = variants.back();
			variants.pop_back(); variants.pop_back();
			variants.push_back(varB_copy);
		    }
		}
	    }
        }
    }

    return variantAlleles;
}

ostream& operator<<(ostream& out, VariantAllele& var) {
    out << var.position << " " << var.ref << " -> " << var.alt;
    return out;
}

bool operator<(const VariantAllele& a, const VariantAllele& b) {
    return a.repr < b.repr;
}

map<pair<int, int>, int> Variant::getGenotypeIndexesDiploid(void) {

    map<pair<int, int>, int> genotypeIndexes;
    //map<int, map<Genotype*, int> > vcfGenotypeOrder;
    vector<int> indexes;
    for (int i = 0; i < alleles.size(); ++i) {
        indexes.push_back(i);
    }
    int ploidy = 2; // ONLY diploid
    vector<vector<int> > genotypes = multichoose(ploidy, indexes);
    for (vector<vector<int> >::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        sort(g->begin(), g->end());  // enforce e.g. 0/1, 0/2, 1/2 ordering over reverse
        // XXX this does not handle non-diploid!!!!
        int j = g->front();
        int k = g->back();
        genotypeIndexes[make_pair(j, k)] = (k * (k + 1) / 2) + j;
    }
    return genotypeIndexes;

}

void Variant::updateAlleleIndexes(void) {
    // adjust the allele index
    altAlleleIndexes.clear();
    int m = 0;
    for (vector<string>::iterator a = alt.begin();
            a != alt.end(); ++a, ++m) {
        altAlleleIndexes[*a] = m;
    }
}

// TODO only works on "A"llele variant fields
void Variant::removeAlt(string& altAllele) {

    int altIndex = getAltAlleleIndex(altAllele);  // this is the alt-relative index, 0-based

    for (map<string, int>::iterator c = vcf->infoCounts.begin(); c != vcf->infoCounts.end(); ++c) {
        int count = c->second;
        if (count == ALLELE_NUMBER) {
            string key = c->first;
            map<string, vector<string> >::iterator v = info.find(key);
            if (v != info.end()) {
                vector<string>& vals = v->second;
                vector<string> tokeep;
                int i = 0;
                for (vector<string>::iterator a = vals.begin(); a != vals.end(); ++a, ++i) {
                    if (i != altIndex) {
                        tokeep.push_back(*a);
                    }
                }
                vals = tokeep;
            }
        }
    }

    for (map<string, int>::iterator c = vcf->formatCounts.begin(); c != vcf->formatCounts.end(); ++c) {
        int count = c->second;
        if (count == ALLELE_NUMBER) {
            string key = c->first;
	    for (map<string, map<string, vector<string> > >::iterator s = samples.begin(); s != samples.end(); ++s) {
		map<string, vector<string> >& sample = s->second;
		map<string, vector<string> >::iterator v = sample.find(key);
		if (v != sample.end()) {
		    vector<string>& vals = v->second;
		    vector<string> tokeep;
		    int i = 0;
		    for (vector<string>::iterator a = vals.begin(); a != vals.end(); ++a, ++i) {
			if (i != altIndex) {
			    tokeep.push_back(*a);
			}
		    }
		    vals = tokeep;
		}
	    }
        }
    }

    int altSpecIndex = altIndex + 1; // this is the genotype-spec index, ref=0, 1-based for alts

    vector<string> newalt;
    map<int, int> alleleIndexMapping;
    // setup the new alt string
    alleleIndexMapping[0] = 0; // reference allele remains the same
    int i = 1; // current index
    int j = 1; // new index
    for (vector<string>::iterator a = alt.begin(); a != alt.end(); ++a, ++i) {
        if (i != altSpecIndex) {
            newalt.push_back(*a);
            // get the mapping between new and old allele indexes
            alleleIndexMapping[i] = j;
            ++j;
        } else {
            alleleIndexMapping[i] = NULL_ALLELE;
        }
    }

    // fix the sample genotypes, removing reference to the old allele
    for (map<string, map<string, vector<string> > >::iterator s = samples.begin(); s != samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<int, int> genotype = decomposeGenotype(sample["GT"].front());
        map<int, int> newGenotype;
        for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
            newGenotype[alleleIndexMapping[g->first]] += g->second;
        }
        sample["GT"].clear();
        sample["GT"].push_back(genotypeToString(newGenotype));
    }

    // reset the alt
    alt = newalt;

    // and the alleles
    alleles.clear();
    alleles.push_back(ref);
    alleles.insert(alleles.end(), alt.begin(), alt.end());

    updateAlleleIndexes();

}

// union of lines in headers of input files
string unionInfoHeaderLines(string& s1, string& s2) {
    vector<string> lines1 = split(s1, "\n");
    vector<string> lines2 = split(s2, "\n");
    vector<string> result;
    set<string> l2;
    string lastHeaderLine; // this one needs to be at the end
    for (vector<string>::iterator s = lines2.begin(); s != lines2.end(); ++s) {
	if (s->substr(0,6) == "##INFO") {
	    l2.insert(*s);
	}
    }
    for (vector<string>::iterator s = lines1.begin(); s != lines1.end(); ++s) {
	if (l2.count(*s)) {
	    l2.erase(*s);
	}
	if (s->substr(0,6) == "#CHROM") {
	    lastHeaderLine = *s;
	} else {
	    result.push_back(*s);
	}
    }
    for (set<string>::iterator s = l2.begin(); s != l2.end(); ++s) {
	result.push_back(*s);
    }
    if (lastHeaderLine.empty()) {
	cerr << "could not find CHROM POS ... header line" << endl;
	exit(1);
    }
    result.push_back(lastHeaderLine);
    return join(result, "\n");
}

} // end namespace vcf
