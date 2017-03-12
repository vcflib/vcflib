#include "Variant.h"
#include <utility>


namespace vcflib {



void Variant::parse(string& line, bool parseSamples) {

    // clean up potentially variable data structures
    info.clear();
    infoFlags.clear();
    format.clear();
    alt.clear();
    alleles.clear();

    // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT [SAMPLE1 .. SAMPLEN]
    vector<string> fields = split(line, '\t');
    if (fields.size() < 7) {
        cerr << "broken VCF record (less than 7 fields)" << endl
             << line << endl;
        exit(1);
    }

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
    int n = 0;
    for (vector<string>::iterator a = alt.begin();
            a != alt.end(); ++a, ++n) {
        altAlleleIndexes[*a] = n;
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
                split(kv.at(1), ',', info[kv.at(0)]);
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
        for (; sample != fields.end() && sampleName != sampleNames.end(); 
	     ++sample, ++sampleName) {
	  string& name = *sampleName;
	  
	  vector<string> samplefields = split(*sample, ':');
	  vector<string>::iterator i = samplefields.begin();
	  
	  for (vector<string>::iterator f = format.begin(); 
	       f != format.end(); ++f) {
	    
	    if(i != samplefields.end()){
	      samples[name][*f] = split(*i, ','); ++i;
	    }
	    else{
	      std::vector<string> missing;
	      missing.push_back(".");
	      samples[name][*f] = missing;
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
            cerr << *sample << endl;
            exit(1);
        }
    } else if (!parseSamples) {
        originalLine = line;
    }

    //return true; // we should be catching exceptions...
}

bool Variant::is_sv(){
    return this->info.find("SVTYPE") != this->info.end();
}

bool Variant::canonicalize_sv(FastaReference& fasta_reference, vector<FastaReference*> insertions, int max_interval){
    
            bool variant_acceptable = true;
            bool do_external_insertions = !insertions.empty();
            int32_t sv_len = 0;
            bool var_is_sv = false;
            FastaReference* insertion_fasta;

#ifdef DEBUG
            if (do_external_insertions){
                insertion_fasta = insertions[0];
                for (auto x : insertion_fasta->index->sequenceNames){
                    cerr << x << endl;
                }
            }
#endif


    std::function<bool(const string&)> allATGC = [](const string& s){
    for (string::const_iterator c = s.begin(); c != s.end(); ++c) {
        char b = *c;
        if (b != 'A' && b != 'T' && b != 'G' && b != 'C') {
            return false;
        }
    }
    return true;
    };

        if (!this->is_sv()){
            return true;
        }


                if (this->info.find("SVTYPE") != this->info.end() &&
                        !(this->info["SVTYPE"].empty())){

                    var_is_sv = true;

                    for (int alt_pos = 0; alt_pos < this->alt.size(); ++alt_pos) {
                        string a = this->alt[alt_pos];
                            // These should be normalized-ish
                            // Ref field might be "N"
                            // Alt field could be <INS>, but it *should* be the inserted sequence,
                            // but that might be tucked away in an info field
                            //
                            //
                            // From the VCF4.2 Specs:
                            // END is POS + LEN(REF) - 1
                            // SVLEN is the difference in length between REF and ALT alleles.


                        if (!(this->info["SVTYPE"][alt_pos] == "INV" ||
                                        this->info["SVTYPE"][alt_pos] == "DEL" ||
                                        this->info["SVTYPE"][alt_pos] == "INS") || this->alt.size() > 1){
                                variant_acceptable = false;
                                break;                                                                                                                                  
                        }

                        if (this->info.find("SVLEN") != this->info.end()){
                            int32_t pre_abs = stoi(this->info["SVLEN"][alt_pos]); 
                            sv_len = abs(pre_abs);
                        }
                        else if (this->info.find("END") != this->info.end()){
                            int32_t pre_abs = stoi(this->info["END"][alt_pos]) - (this->position);
                            sv_len = abs( pre_abs );

                        }
                        else{
                            // If we have neither, we'll ignore it.
                            variant_acceptable = false;
                            break;
                        }

                        if (this->info["SVTYPE"][alt_pos] == "INS" || a == "<INS>"){
                            this->ref.assign(fasta_reference.getSubSequence(this->sequenceName, this->position, 1));
                            //if (this->alt[alt_pos] == "<INS>"){
                            //    variant_acceptable = false;
                            //}
                            if (do_external_insertions){
                                insertion_fasta = insertions[0];
                                string var_name;
                                if (alt[alt_pos][0] == '<' && alt[alt_pos][ alt[alt_pos].size() - 1 ] == '>'){
                                    string shortname (alt[alt_pos], 1, alt[alt_pos].length() - 2 );
                                    var_name = shortname;
                                }
                                else{
                                    var_name = alt[alt_pos]; 
                                }
                                #ifdef DEBUG
                                    cerr << "My variant name is: " << var_name << endl;
                                #endif
                                if (insertion_fasta->index->find(var_name) != insertion_fasta->index->end()){
                                    this->ref.assign(fasta_reference.getSubSequence(this->sequenceName, this->position, 1 ));
                                    #ifdef DEBUG
                                        cerr << "Replacing insertion with sequence of " << var_name << endl;
                                    #endif
                                    this->alt[alt_pos] = fasta_reference.getSubSequence(this->sequenceName, this->position, 1) + insertion_fasta->getSequence(var_name);
                                    this->updateAlleleIndexes();
                                }
                            }
                            else if (allATGC(a)){
                                // we don't have to do anything - our INS is already in the correct format!
                            }
                            else{
                                variant_acceptable = false;
                            }
                        }
                        else if (a == "<DEL>" || this->info["SVTYPE"][alt_pos] == "DEL"){


                            this->ref.assign(fasta_reference.getSubSequence(this->sequenceName, this->position, sv_len + 1 ));

                            this->alt[alt_pos].assign(fasta_reference.getSubSequence(this->sequenceName, this->position, 1));

                            if (this->ref.size() != sv_len + 1){
                                cerr << "Variant made is incorrect size" << endl;
                                cerr << this->ref.size() - 1 << "\t" << sv_len << endl;
                                cerr << this->ref[this->ref.size() - 1] << "\t" << endl;
                                variant_acceptable = false;
                            }
                            this->updateAlleleIndexes();

                        }
                        else if (a == "<INV>" || this->info["SVTYPE"][alt_pos] == "INV"){
                            this->ref = fasta_reference.getSubSequence(this->sequenceName, this->position, sv_len);
                            string alt_str(fasta_reference.getSubSequence(this->sequenceName, this->position, sv_len));
                            reverse(alt_str.begin(), alt_str.end());
                            this->alt[alt_pos] = alt_str;

                            // add 3 bases padding to right side 
                            //this->ref.insert(0, reference.getSubSequence(this->sequenceName, this->position - 3, 3));
                            //this->alt[alt_pos].insert(0, reference.getSubSequence(this->sequenceName, this->position - 3, 3));
                            //this->position = this->position - 3;
                            this->updateAlleleIndexes();

                            //variant_acceptable = false;

                        }
                        else{
                            variant_acceptable = false;
                        }
                    }

                }


    return variant_acceptable;
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
        if (key == "FILTER") { // hack to use FILTER as an "info" field (why the hack?)
            return FIELD_STRING;
        }
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
            if (key == "FILTER") {
              return filter;
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
                exit(1);
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
        // ensure there are no empty fields
        if (var.sequenceName.empty()) var.sequenceName = ".";
        if (var.id.empty()) var.id = ".";
        if (var.ref.empty()) var.ref = ".";
        if (var.alt.empty()) var.alt.push_back(".");
        if (var.filter.empty()) var.filter = ".";

        out << var.sequenceName << "\t"
            << var.position << "\t"
            << var.id << "\t"
            << var.ref << "\t";
        // report the list of alternate alleles.
        var.printAlt(out);
        out << "\t"
            << var.quality << "\t"
            << var.filter << "\t";
        if (var.info.empty() && var.infoFlags.empty()) {
            out << ".";
        } else {
            for (map<string, vector<string> >::iterator i = var.info.begin(); i != var.info.end(); ++i) {
                if (!i->second.empty()) {
                    out << ((i == var.info.begin()) ? "" : ";") << i->first << "=" << join(i->second, ",");
                }
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
                            if (g != sample.end() && !g->second.empty()) {
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
                } else if (tokenstr == "FILTER") {
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
        for (unsigned int i = 0; i <  filterspec.size(); ++i) {
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
            RuleToken token = rulesCopy.front();
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
                VariantFieldType vtype;
                if (sample.empty()) { // means we are record-specific
                    vtype = var.infoType(token.value);
                } else {
                    vtype = var.formatType(token.value);
                    //cout << "type = " << type << endl;
                }
                //cout << "type: " << type << endl;

                if (vtype == FIELD_INTEGER || vtype == FIELD_FLOAT) {
                    token.type = RuleToken::NUMERIC_VARIABLE;
                    token.number = var.getValueFloat(token.value, sample, index);
                    //cerr << "number: " << token.number << endl;
                } else if (vtype == FIELD_BOOL) {
                    token.type = RuleToken::BOOLEAN_VARIABLE;
                    token.state = var.getValueBool(token.value, sample, index);
                    //cerr << "state: " << token.state << endl;
                } else if (vtype == FIELD_STRING) {
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
                    } else if (a.type == RuleToken::STRING_VARIABLE && b.type == RuleToken::NUMERIC_VARIABLE) {
                        r.state = (convert(b.number) == a.str);
                    } else if (b.type == RuleToken::STRING_VARIABLE && a.type == RuleToken::NUMERIC_VARIABLE) {
                        r.state = (convert(a.number) == b.str);
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

void VariantFilter::removeFilteredGenotypes(Variant& var, bool keepInfo) {

    for (vector<string>::iterator s = var.sampleNames.begin(); s != var.sampleNames.end(); ++s) {
        string& name = *s;
        if (!passes(var, name)) {
        	if (keepInfo) {
				var.samples[name]["GT"].clear();
				var.samples[name]["GT"].push_back("./.");
        	}
        	else {
			    var.samples.erase(name);
        	}
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

// non-destructive version of above
string VariantCallFile::headerWithSampleNames(vector<string>& newSamples) {
    // regenerate the last line of the header
    if (newSamples.empty()) return header;
    vector<string> headerLines = split(header, '\n');
    vector<string> colnames = split(headerLines.at(headerLines.size() - 1), '\t'); // get the last, update the samples
    vector<string> newcolnames;
    unsigned int colCount = colnames.size(); // used to be hard-coded 9, hopefully the dynamic colCount isn't an issue
    if (colCount < 8)
    {
        cout << "VCF file is not suitable for use because it does not have a format field." << endl;
        exit(0);
    }
    newcolnames.resize(colCount + newSamples.size());
    copy(colnames.begin(), colnames.begin() + colCount, newcolnames.begin());
    copy(newSamples.begin(), newSamples.end(), newcolnames.begin() + colCount);
    headerLines.at(headerLines.size() - 1) = join(newcolnames, "\t");
    return join(headerLines, "\n");
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
    string id = "ID=" + tag + ",";
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
    string id = "ID=" + tag + ",";
    for (vector<string>::iterator s = headerLines.begin(); s != headerLines.end(); ++s) {
        string& headerLine = *s;
        if (headerLine.find("##FORMAT") == 0) {
            if (headerLine.find(id) == string::npos) {
                newHeader.push_back(headerLine);
            }
        } else {
            newHeader.push_back(headerLine);
        }
    }
    header = join(newHeader, "\n");
}

vector<string> VariantCallFile::getHeaderLinesFromFile()
{
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
                    return vector<string>();
                }
                firstRecord = true;
                break;
            }
        }
    }
    return split(headerStr, "\n");
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
                    return false;
                }
                firstRecord = true;
                break;
            }
        }
    }
    this->vcf_header = headerStr;

    return parseHeader(headerStr);

}

bool VariantCallFile::parseHeader(string& hs) {

    if (hs.empty()) return false;
    if (hs.substr(hs.size() - 1, 1) == "\n") {
	hs.erase(hs.size() - 1, 1); // remove trailing newline
    }
    header = hs; // stores the header in the object instance

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
        if (firstRecord && !justSetRegion) {
            if (!line.empty() && line.substr(0,1) != "#") {
                var.parse(line, parseSamples);
                firstRecord = false;
                _done = false;
                return true;
            } else {
                return false;
            }
        }
        if (usingTabix) {
            if (justSetRegion && !line.empty() && line.substr(0,1) != "#") {
                if (firstRecord) {
                    firstRecord = false;
                }
                var.parse(line, parseSamples);
                line.clear();
                justSetRegion = false;
                _done = false;
                return true;
            } else if (tabixFile->getNextLine(line)) {
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
    return setRegion(regionstr.str());
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
	    justSetRegion = true;
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

map<int, int> decomposeGenotype(const string& genotype) {
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

vector<int> decomposePhasedGenotype(const string& genotype) {
    string splitter = "/";
    if (genotype.find("|") != string::npos) {
        splitter = "|";
    }
    vector<string> haps = split(genotype, splitter);
    if (haps.size() > 1 && splitter == "/") {
        cerr << "could not find '|' in genotype, cannot decomposePhasedGenotype on unphased genotypes" << endl;
        exit(1);
    }
    vector<int> decomposed;
    for (vector<string>::iterator h = haps.begin(); h != haps.end(); ++h) {
        int alt;
        if (*h == ".") {
            decomposed.push_back(NULL_ALLELE);
        } else {
            convert(*h, alt);
            decomposed.push_back(alt);
        }
    }
    return decomposed;
}

string genotypeToString(const map<int, int>& genotype) {
    vector<int> s;
    for (map<int, int>::const_iterator g = genotype.begin(); g != genotype.end(); ++g) {
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

string phasedGenotypeToString(const vector<int>& genotype) {
    vector<string> r;
    for (vector<int>::const_iterator i = genotype.begin(); i != genotype.end(); ++i) {
        if (*i == NULL_ALLELE) r.push_back(".");
        else r.push_back(convert(*i));
    }
    return join(r, "|");
}

bool isHet(const map<int, int>& genotype) {
    return genotype.size() > 1;
}

bool isHom(const map<int, int>& genotype) {
    return genotype.size() == 1;
}

bool hasNonRef(const map<int, int>& genotype) {
    for (map<int, int>::const_iterator g = genotype.begin(); g != genotype.end(); ++g) {
        if (g->first != 0) {
            return true;
        }
    }
    return false;
}

bool isHomRef(const map<int, int>& genotype) {
    return isHom(genotype) && !hasNonRef(genotype);
}

bool isHomNonRef(const map<int, int>& genotype) {
    return isHom(genotype) && hasNonRef(genotype);
}

bool isNull(const map<int, int>& genotype) {
    return genotype.find(NULL_ALLELE) != genotype.end();
}

int ploidy(const map<int, int>& genotype) {
    int i = 0;
    for (map<int, int>::const_iterator g = genotype.begin(); g != genotype.end(); ++g) {
        i += g->second;
    }
    return i;
}

// generates cigar from allele parsed by parsedAlternates
string varCigar(vector<VariantAllele>& vav, bool xForMismatch) {
    string cigar;
    pair<int, string> element;
    for (vector<VariantAllele>::iterator v = vav.begin(); v != vav.end(); ++v) {
        VariantAllele& va = *v;
        if (va.ref != va.alt) {
            if (element.second == "M") {
                cigar += convert(element.first) + element.second;
                element.second = ""; element.first = 0;
            }
            if (va.ref.size() == va.alt.size()) {
                cigar += convert(va.ref.size()) + (xForMismatch ? "X" : "M");
            } else if (va.ref.size() > va.alt.size()) {
                cigar += convert(va.ref.size() - va.alt.size()) + "D";
            } else {
                cigar += convert(va.alt.size() - va.ref.size()) + "I";
            }
        } else {
            if (element.second == "M") {
                element.first += va.ref.size();
            } else {
                element = make_pair(va.ref.size(), "M");
            }
        }
    }
    if (element.second == "M") {
        cigar += convert(element.first) + element.second;
    }
    element.second = ""; element.first = 0;
    return cigar;
}

map<string, vector<VariantAllele> > Variant::parsedAlternates(bool includePreviousBaseForIndels,
                                                              bool useMNPs,
                                                              bool useEntropy,
                                                              float matchScore,
                                                              float mismatchScore,
                                                              float gapOpenPenalty,
                                                              float gapExtendPenalty,
                                                              float repeatGapExtendPenalty,
                                                              string flankingRefLeft,
                                                              string flankingRefRight) {

    map<string, vector<VariantAllele> > variantAlleles;

    // add the reference allele
    variantAlleles[ref].push_back(VariantAllele(ref, ref, position));

    // single SNP case, no ambiguity possible, no need to spend a lot of
    // compute aligning ref and alt fields
    if (alt.size() == 1 && ref.size() == 1 && alt.front().size() == 1) {
        variantAlleles[alt.front()].push_back(VariantAllele(ref, alt.front(), position));
        return variantAlleles;
    }

    // padding is used to ensure a stable alignment of the alternates to the reference
    // without having to go back and look at the full reference sequence
    int paddingLen = max(10, (int) (ref.size()));  // dynamically determine optimum padding length
    for (vector<string>::iterator a = alt.begin(); a != alt.end(); ++a) {
        string& alternate = *a;
        paddingLen = max(paddingLen, (int) (alternate.size()));
    }
    char padChar = 'Z';
    char anchorChar = 'Q';
    string padding(paddingLen, padChar);

    // this 'anchored' string is done for stability
    // the assumption is that there should be a positional match in the first base
    // this is true for VCF 4.1, and standard best practices
    // using the anchor char ensures this without other kinds of realignment
    string reference_M;
    if (flankingRefLeft.empty() && flankingRefRight.empty()) {
        reference_M = padding + ref + padding;
        reference_M[paddingLen] = anchorChar;
    } else {
        reference_M = flankingRefLeft + ref + flankingRefRight;
        paddingLen = flankingRefLeft.size();
    }

    // passed to sw.Align
    unsigned int referencePos;

    string cigar;
    
    for (vector<string>::iterator a = alt.begin(); a != alt.end(); ++a) {
      
      string& alternate = *a;
      vector<VariantAllele>& variants = variantAlleles[alternate];
      string alternateQuery_M;
      if (flankingRefLeft.empty() && flankingRefRight.empty()) {
	alternateQuery_M = padding + alternate + padding;
	alternateQuery_M[paddingLen] = anchorChar;
      } else {
	alternateQuery_M = flankingRefLeft + alternate + flankingRefRight;
      }
      //const unsigned int alternateLen = alternate.size();
      
      if (true) {
	CSmithWatermanGotoh sw(matchScore, 
			       mismatchScore, 
			       gapOpenPenalty, 
			       gapExtendPenalty);
	if (useEntropy) sw.EnableEntropyGapPenalty(1);
	if (repeatGapExtendPenalty != 0){
	  sw.EnableRepeatGapExtensionPenalty(repeatGapExtendPenalty);
	}
	sw.Align(referencePos, cigar, reference_M, alternateQuery_M);
      } else {  // disabled for now
	StripedSmithWaterman::Aligner aligner;
	StripedSmithWaterman::Filter sswFilter;
	StripedSmithWaterman::Alignment alignment;
	aligner.Align(alternateQuery_M.c_str(), 
		      reference_M.c_str(), 
		      reference_M.size(), sswFilter, &alignment);
	cigar = alignment.cigar_string;
      }

      // left-realign the alignment...
      
      vector<pair<int, string> > cigarData = splitCigar(cigar);
      
      if (cigarData.front().second != "M" 
	  || cigarData.back().second != "M"
	  || cigarData.front().first < paddingLen 
	  || cigarData.back().first < paddingLen) {
	cerr << "parsedAlternates: alignment does not start with match over padded sequence" << endl;
	cerr << cigar << endl;
	cerr << reference_M << endl;
	cerr << alternateQuery_M << endl;
	exit(1);
      } else {
	cigarData.front().first -= paddingLen;
	cigarData.back().first -= paddingLen;;
      }
      //cigarData = cleanCigar(cigarData);
      cigar = joinCigar(cigarData);
      
      int altpos = 0;
      int refpos = 0;
      
      for (vector<pair<int, string> >::iterator e = cigarData.begin(); 
	   e != cigarData.end(); ++e) {
	
	int len = e->first;
	string type = e->second;
	
	switch (type.at(0)) {
	case 'I':
	  if (includePreviousBaseForIndels) {
	    if (!variants.empty() && 
		variants.back().ref != variants.back().alt) {
	      VariantAllele a = 
		VariantAllele("", 
			      alternate.substr(altpos, len), 
			      refpos + position);
	      variants.back() = variants.back() + a;
	    } else {
	      VariantAllele a = 
		VariantAllele(ref.substr(refpos - 1, 1),
			      alternate.substr(altpos - 1, len + 1),
			      refpos + position - 1);
	      variants.push_back(a);
	    }
	  } else {
	    variants.push_back(VariantAllele("", 
					     alternate.substr(altpos, len),
					     refpos + position));
	  }
	  altpos += len;
	  break;
	case 'D':
	  if (includePreviousBaseForIndels) {
	    if (!variants.empty() &&
		variants.back().ref != variants.back().alt) {
	      VariantAllele a 
		= VariantAllele(ref.substr(refpos, len)
				, "", refpos + position);
	      variants.back() = variants.back() + a;
	      } else {
	      VariantAllele a 
		= VariantAllele(ref.substr(refpos - 1, len + 1),
				alternate.substr(altpos - 1, 1),
				refpos + position - 1);
	      variants.push_back(a);
	    }
	  } else {
	    variants.push_back(VariantAllele(ref.substr(refpos, len), 
					     "", refpos + position));
	  }
	  refpos += len;
	  break;

	  // zk has added (!variants.empty()) solves the seg fault in 
          // vcfstats, but need to test
	case 'M':
	  {
	    for (int i = 0; i < len; ++i) {
	      VariantAllele a 
		= VariantAllele(ref.substr(refpos + i, 1),
				alternate.substr(altpos + i, 1),
				(refpos + i + position));
	      if (useMNPs && (!variants.empty()) &&
		  variants.back().ref.size() == variants.back().alt.size()
		  && variants.back().ref != variants.back().alt) {
		  variants.back() = variants.back() + a;
	      } else {
		variants.push_back(a);
	      }
	    }
	  }
	  refpos += len;
	  altpos += len;
	  break;
	case 'S':
	  {
	    refpos += len;
	    altpos += len;
	    break;
	  }
	default:
	  {
	    break;
	  }
	}
      }
    }
    return variantAlleles;
}

map<string, vector<VariantAllele> > Variant::flatAlternates(void) {
    map<string, vector<VariantAllele> > variantAlleles;
    for (vector<string>::iterator a = alt.begin(); a != alt.end(); ++a) {
        string& alternate = *a;
        vector<VariantAllele>& variants = variantAlleles[alternate];
        variants.push_back(VariantAllele(ref, alternate, position));
    }
    return variantAlleles;
}

set<string> Variant::altSet(void) {
    set<string> altset(alt.begin(), alt.end());
    return altset;
}

ostream& operator<<(ostream& out, VariantAllele& var) {
    out << var.position << " " << var.ref << " -> " << var.alt;
    return out;
}

VariantAllele operator+(const VariantAllele& a, const VariantAllele& b) {
    return VariantAllele(a.ref + b.ref, a.alt + b.alt, a.position);
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
    
    for (map<string, int>::iterator c = vcf->infoCounts.begin(); 
	 c != vcf->infoCounts.end(); ++c) {
      int count = c->second;
      if (count == ALLELE_NUMBER) {
	string key = c->first;
	map<string, vector<string> >::iterator v = info.find(key);
	if (v != info.end()) {
	  vector<string>& vals = v->second;
	  vector<string> tokeep;
	  int i = 0;
	  for (vector<string>::iterator a = vals.begin(); 
	       a != vals.end(); ++a, ++i) {
	    if (i != altIndex) {
	      tokeep.push_back(*a);
	    }
	  }
	  vals = tokeep;
	}
      }
    }
    
    for (map<string, int>::iterator c = vcf->formatCounts.begin(); 
	 c != vcf->formatCounts.end(); ++c) {
      int count = c->second;
      if (count == ALLELE_NUMBER) {
            string key = c->first;
            for (map<string, map<string, vector<string> > >::iterator 
		   s = samples.begin(); s != samples.end(); ++s) {
	      map<string, vector<string> >& sample = s->second;
	      map<string, vector<string> >::iterator v = sample.find(key);
	      if (v != sample.end()) {
		vector<string>& vals = v->second;
		vector<string> tokeep;
		int i = 0;
		for (vector<string>::iterator a = vals.begin(); 
		     a != vals.end(); ++a, ++i) {
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
    alleleIndexMapping[NULL_ALLELE] = NULL_ALLELE; // null allele remains the same
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
    map<string, int> samplePloidy;
    for (map<string, map<string, vector<string> > >::iterator s = samples.begin(); s != samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        if (sample.find("GT") != sample.end()) {
            string& gt = sample["GT"].front();
            string splitter = "/";
            if (gt.find("|") != string::npos) {
                splitter = "|";
            }

            if (splitter == "/") {
                samplePloidy[s->first] = split(gt, splitter).size();
                map<int, int> genotype = decomposeGenotype(sample["GT"].front());
                map<int, int> newGenotype;
                for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
                    newGenotype[alleleIndexMapping[g->first]] += g->second;
                }
                sample["GT"].clear();
                sample["GT"].push_back(genotypeToString(newGenotype));
            } else {
                samplePloidy[s->first] = split(gt, splitter).size();
                vector<int> genotype = decomposePhasedGenotype(sample["GT"].front());
                vector<int> newGenotype;
                for (vector<int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
                    newGenotype.push_back(alleleIndexMapping[*g]);
                }
                sample["GT"].clear();
                sample["GT"].push_back(phasedGenotypeToString(newGenotype));
            }
        }
    }

    set<int> ploidies;
    for (map<string, int>::iterator p = samplePloidy.begin(); p != samplePloidy.end(); ++p) {
        ploidies.insert(p->second);
    }

    // fix the sample genotype likelihoods, removing reference to the old allele
    // which GL fields should we remove?
    vector<int> toRemove;
    toRemove.push_back(altSpecIndex);
    map<int, map<int, int> > glMappingByPloidy;
    for (set<int>::iterator p = ploidies.begin(); p != ploidies.end(); ++p) {
        glMappingByPloidy[*p] = glReorder(*p, alt.size() + 1, alleleIndexMapping, toRemove);
    }

    for (map<string, map<string, vector<string> > >::iterator s = samples.begin(); s != samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator glsit = sample.find("GL");
        if (glsit != sample.end()) {
            vector<string>& gls = glsit->second; // should be split already
            map<int, string> newgls;
            map<int, int>& newOrder = glMappingByPloidy[samplePloidy[s->first]];
            int i = 0;
            for (vector<string>::iterator g = gls.begin(); g != gls.end(); ++g, ++i) {
                int j = newOrder[i];
                if (j != -1) {
                    newgls[i] = *g;
                }
            }
            // update the gls
            gls.clear();
            for (map<int, string>::iterator g = newgls.begin(); g != newgls.end(); ++g) {
                gls.push_back(g->second);
            }
        }
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

string mergeCigar(const string& c1, const string& c2) {
    vector<pair<int, string> > cigar1 = splitCigar(c1);
    vector<pair<int, string> > cigar2 = splitCigar(c2);
    // check if the middle elements are the same
    if (cigar1.back().second == cigar2.front().second) {
        cigar1.back().first += cigar2.front().first;
        cigar2.erase(cigar2.begin());
    }
    for (vector<pair<int, string> >::iterator c = cigar2.begin(); c != cigar2.end(); ++c) {
        cigar1.push_back(*c);
    }
    return joinCigar(cigar1);
}

vector<pair<int, string> > splitCigar(const string& cigarStr) {
    vector<pair<int, string> > cigar;
    string number;
    string type;
    // strings go [Number][Type] ...
    for (string::const_iterator s = cigarStr.begin(); s != cigarStr.end(); ++s) {
        char c = *s;
        if (isdigit(c)) {
            if (type.empty()) {
                number += c;
            } else {
                // signal for next token, push back the last pair, clean up
                cigar.push_back(make_pair(atoi(number.c_str()), type));
                number.clear();
                type.clear();
                number += c;
            }
        } else {
            type += c;
        }
    }
    if (!number.empty() && !type.empty()) {
        cigar.push_back(make_pair(atoi(number.c_str()), type));
    }
    return cigar;
}

list<pair<int, string> > splitCigarList(const string& cigarStr) {
    list<pair<int, string> > cigar;
    string number;
    string type;
    // strings go [Number][Type] ...
    for (string::const_iterator s = cigarStr.begin(); s != cigarStr.end(); ++s) {
        char c = *s;
        if (isdigit(c)) {
            if (type.empty()) {
                number += c;
            } else {
                // signal for next token, push back the last pair, clean up
                cigar.push_back(make_pair(atoi(number.c_str()), type));
                number.clear();
                type.clear();
                number += c;
            }
        } else {
            type += c;
        }
    }
    if (!number.empty() && !type.empty()) {
        cigar.push_back(make_pair(atoi(number.c_str()), type));
    }
    return cigar;
}

vector<pair<int, string> > cleanCigar(const vector<pair<int, string> >& cigar) {
    vector<pair<int, string> > cigarClean;
    for (vector<pair<int, string> >::const_iterator c = cigar.begin(); c != cigar.end(); ++c) {
        if (c->first > 0) {
            cigarClean.push_back(*c);
        }
    }
    return cigarClean;
}

string joinCigar(const vector<pair<int, string> >& cigar) {
    string cigarStr;
    for (vector<pair<int, string> >::const_iterator c = cigar.begin(); c != cigar.end(); ++c) {
        if (c->first) {
            cigarStr += convert(c->first) + c->second;
        }
    }
    return cigarStr;
}

string joinCigar(const vector<pair<int, char> >& cigar) {
    string cigarStr;
    for (vector<pair<int, char> >::const_iterator c = cigar.begin(); c != cigar.end(); ++c) {
        if (c->first) {
            cigarStr += convert(c->first) + string(1, c->second);
        }
    }
    return cigarStr;
}

string joinCigarList(const list<pair<int, string> >& cigar) {
    string cigarStr;
    for (list<pair<int, string> >::const_iterator c = cigar.begin(); c != cigar.end(); ++c) {
        cigarStr += convert(c->first) + c->second;
    }
    return cigarStr;
}

int cigarRefLen(const vector<pair<int, char> >& cigar) {
    int len = 0;
    for (vector<pair<int, char> >::const_iterator c = cigar.begin(); c != cigar.end(); ++c) {
        if (c->second == 'M' || c->second == 'D' || c->second == 'X') {
            len += c->first;
        }
    }
    return len;
}

int cigarRefLen(const vector<pair<int, string> >& cigar) {
    int len = 0;
    for (vector<pair<int, string> >::const_iterator c = cigar.begin(); c != cigar.end(); ++c) {
        if (c->second == "M" || c->second == "D" || c->second == "X") {
            len += c->first;
        }
    }
    return len;
}

bool isEmptyCigarElement(const pair<int, string>& elem) {
    return elem.first == 0;
}

list<list<int> > _glorder(int ploidy, int alts) {
    if (ploidy == 1) {
        list<list<int> > results;
        for (int n = 0; n < alts; ++n) {
            list<int> v;
            v.push_back(n);
            results.push_back(v);
        }
        return results;
    } else {
        list<list<int> > results;
        for (int n = 0; n < alts; ++n) {
            list<list<int> > x = _glorder(ploidy - 1, alts);
            for (list<list<int> >::iterator v = x.begin(); v != x.end(); ++v) {
                if (v->front() <= n) {
                    v->push_front(n);
                    results.push_back(*v);
                }
            }
        }
        return results;
    }
}

// genotype likelihood-ordering of genotypes, where each genotype is a
// list of integers (as written in the GT field)
list<list<int> > glorder(int ploidy, int alts) {
    list<list<int> > results = _glorder(ploidy, alts);
    for (list<list<int> >::iterator v = results.begin(); v != results.end(); ++v) {
        v->reverse();
    }
    return results;
}

// which genotype likelihoods would include this alternate allele
list<int> glsWithAlt(int alt, int ploidy, int numalts) {
    list<int> gls;
    list<list<int> > orderedGenotypes = glorder(ploidy, numalts);
    int i = 0;
    for (list<list<int> >::iterator v = orderedGenotypes.begin(); v != orderedGenotypes.end(); ++v, ++i) {
        for (list<int>::iterator q = v->begin(); q != v->end(); ++q) {
            if (*q == alt) {
                gls.push_back(i);
                break;
            }
        }
    }
    return gls;
}

// describes the mapping between the old gl ordering and and a new
// one in which the GLs including the old alt have been removed
// a map to -1 means "remove"
map<int, int> glReorder(int ploidy, int numalts, map<int, int>& alleleIndexMapping, vector<int>& altsToRemove) {
    map<int, int> mapping;
    list<list<int> > orderedGenotypes = glorder(ploidy, numalts);
    for (list<list<int> >::iterator v = orderedGenotypes.begin(); v != orderedGenotypes.end(); ++v) {
        for (list<int>::iterator n = v->begin(); n != v->end(); ++n) {
            *n = alleleIndexMapping[*n];
        }
    }
    list<list<int> > newOrderedGenotypes = glorder(ploidy, numalts - altsToRemove.size());
    map<list<int>, int> newOrderedGenotypesMapping;
    int i = 0;
    // mapping is wrong...
    for (list<list<int> >::iterator v = newOrderedGenotypes.begin(); v != newOrderedGenotypes.end(); ++v, ++i) {
        newOrderedGenotypesMapping[*v] = i;
    }
    i = 0;
    for (list<list<int> >::iterator v = orderedGenotypes.begin(); v != orderedGenotypes.end(); ++v, ++i) {
        map<list<int>, int>::iterator m = newOrderedGenotypesMapping.find(*v);
        if (m != newOrderedGenotypesMapping.end()) {
            //cout << "new gl order of " << i << " is " << m->second << endl;
            mapping[i] = m->second;
        } else {
            //cout << i << " will be removed" << endl;
            mapping[i] = -1;
        }
    }
    return mapping;
}

string Variant::getGenotype(string& sample) {
    map<string, map<string, vector<string> > >::iterator s = samples.find(sample);
    if (s != samples.end()) {
        map<string, vector<string> >::iterator f = s->second.find("GT");
        if (f != s->second.end()) {
            return f->second.front();
        }
    }
    return "";
}

bool Variant::isPhased(void) {
    for (map<string, map<string, vector<string> > >::iterator s = samples.begin(); s != samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator g = sample.find("GT");
        if (g != sample.end()) {
            string gt = g->second.front();
            if (gt.size() > 1 && gt.find("|") == string::npos) {
                return false;
            }
        }
    }
    return true;
}

long Variant::zeroBasedPosition(void) {
    return position - 1;
}

string Variant::vrepr(void) {
    return sequenceName + "\t" + convert(position) + "\t" + join(alleles, ",");
}

// TODO
/*
vector<Variant*> Variant::matchingHaplotypes() {

    int haplotypeStart = var.position;
    int haplotypeEnd = var.position + var.ref.size();

    for (vector<Variant*>::iterator v = overlapping.begin(); v != overlapping.end(); ++v) {
        haplotypeStart = min((*v)->position, (long int) haplotypeStart);
        haplotypeEnd = max((*v)->position + (*v)->ref.size(), (long unsigned int) haplotypeEnd);
    }

    // for everything overlapping and the current variant, construct the local haplotype within the bounds
    // if there is an exact match, the allele in the current VCF does intersect

    string referenceHaplotype = reference.getSubSequence(var.sequenceName, haplotypeStart - 1, haplotypeEnd - haplotypeStart);
    map<string, vector<pair<Variant*, int> > > haplotypes; // map to variant and alt index

    for (vector<Variant*>::iterator v = overlapping.begin(); v != overlapping.end(); ++v) {
        Variant& variant = **v;
        int altindex = 0;
        for (vector<string>::iterator a = variant.alt.begin(); a != variant.alt.end(); ++a, ++altindex) {
            string haplotype = referenceHaplotype;
            // get the relative start and end coordinates for the variant alternate allele
            int relativeStart = variant.position - haplotypeStart;
            haplotype.replace(relativeStart, variant.ref.size(), *a);
            haplotypes[haplotype].push_back(make_pair(*v, altindex));
        }
    }

    Variant originalVar = var;

    // determine the non-intersecting alts
    vector<string> altsToRemove;
    vector<int> altIndexesToRemove;
    for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
        string haplotype = referenceHaplotype;
        int relativeStart = var.position - haplotypeStart;
        haplotype.replace(relativeStart, var.ref.size(), *a);
        map<string, vector<pair<Variant*, int> > >::iterator h = haplotypes.find(haplotype);
        if ((intersecting && !invert && h == haplotypes.end())
            || (intersecting && invert && h != haplotypes.end())
            || (unioning && h != haplotypes.end())) {
            if (tag.empty() && mergeToTag.empty()) {
                altsToRemove.push_back(*a);
            } else {
                if (!tag.empty()) {
                    var.info[tag].push_back(".");
                }
                if (!mergeToTag.empty()) {
                    var.info[mergeToTag].push_back(".");
                }
            }
        } else {
            if (!tag.empty()) {
                var.info[tag].push_back(tagValue);
            }
            // NB: just take the first value for the mergeFromTag
            if (!mergeToTag.empty()) {
                Variant* v = h->second.front().first;
                int index = h->second.front().second;
                if (v->info.find(mergeFromTag) != v->info.end()) {
                    // now you have to find the exact allele...
                    string& otherValue = v->info[mergeFromTag].at(index);
                    var.info[mergeToTag].push_back(otherValue);
                } else if (mergeFromTag == "QUAL") {
                    var.info[mergeToTag].push_back(convert(v->quality));
                } else {
                    var.info[mergeToTag].push_back(".");
                }
            }
        }
    }

    // remove the non-overlapping (intersecting) or overlapping (unioning) alts
    if (intersecting && loci && altsToRemove.size() != var.alt.size()) {
        // we have a match in loci mode, so we should output the whole loci, not just the matching sequence
    } else {
        for (vector<string>::iterator a = altsToRemove.begin(); a != altsToRemove.end(); ++a) {
            var.removeAlt(*a);
        }
    }

    if (unioning) {

        // somehow sort the records and combine them?
        map<long int, vector<Variant*> > variants;
        for (vector<Variant*>::iterator o = overlapping.begin(); o != overlapping.end(); ++o) {
            if ((*o)->position <= var.position && // check ensures proper ordering of variants on output
                outputVariants.find(*o) == outputVariants.end()) {
                outputVariants.insert(*o);
                variants[(*o)->position].push_back(*o);
            }
        }
        // add in the current variant, if it has alts left
        if (!var.alt.empty()) {
            vector<Variant*>& vars = variants[var.position];
            int numalts = 0;
            for (vector<Variant*>::iterator v = vars.begin(); v != vars.end(); ++v) {
                numalts += (*v)->alt.size();
            }
            if (numalts + var.alt.size() == originalVar.alt.size()) {
                variants[var.position].clear();
                variants[var.position].push_back(&originalVar);
            } else {
                variants[var.position].push_back(&var);
            }
        }

        for (map<long int, vector<Variant*> >::iterator v = variants.begin(); v != variants.end(); ++v) {
            for (vector<Variant*>::iterator o = v->second.begin(); o != v->second.end(); ++o) {
                cout << **o << endl;
                lastOutputPosition = max(lastOutputPosition, (*o)->position);
            }
        }
    } else {
        // if any alts remain, output the variant record
        if (!var.alt.empty()) {
            cout << var << endl;
            lastOutputPosition = max(lastOutputPosition, var.position);
        }
    }

}
*/


    VCFHeader::VCFHeader()
    {

        // add manditory fields
        this->header_columns.push_back("#CHROM");
        this->header_columns.push_back("POS");
        this->header_columns.push_back("ID");
        this->header_columns.push_back("REF");
        this->header_columns.push_back("ALT");
        this->header_columns.push_back("QUAL");
        this->header_columns.push_back("FILTER");
        this->header_columns.push_back("INFO");

        // add the line names in order
        // the order is used when outputting as a string
        this->header_line_names_ordered.push_back("##fileFormat");
        this->header_line_names_ordered.push_back("##fileDate");
        this->header_line_names_ordered.push_back("##source");
        this->header_line_names_ordered.push_back("##reference");
        this->header_line_names_ordered.push_back( "##contig");
        this->header_line_names_ordered.push_back("##phasing");
        this->header_line_names_ordered.push_back( "##assembly");

        // add the list names in order
        // the order is used when outputting as a string (getHeaderString)
        this->header_list_names_ordered.push_back("##info");
        this->header_list_names_ordered.push_back("##filter");
        this->header_list_names_ordered.push_back("##format");
        this->header_list_names_ordered.push_back("##alt");
        this->header_list_names_ordered.push_back("##sample");
        this->header_list_names_ordered.push_back("##pedigree");
        this->header_list_names_ordered.push_back("##pedigreedb");

        // initialize the header_lines with the above vector.
        // Set the key as the ##_type_ and the value as an empty string
        // Empty strings are ignored when outputting as string (getHeaderString)
        for (vector<string>::const_iterator header_lines_iter = this->header_line_names_ordered.begin(); header_lines_iter != this->header_line_names_ordered.end(); ++header_lines_iter)
        {
            this->header_lines[(*header_lines_iter)] = "";
        }

        // initialize the header_lines with the above vector.
        // Set the key as the ##_type_ and the value as an empty vector<string>
        // Empty vectors are ignored when outputting as string (getHeaderString)
        for (vector<string>::const_iterator header_lists_iter = this->header_list_names_ordered.begin(); header_lists_iter != this->header_list_names_ordered.end(); ++header_lists_iter)
        {
            this->header_lists[(*header_lists_iter)] = vector<string>(0);
        }

    }

    void VCFHeader::addMetaInformationLine(const string& meta_line)
    {
        // get the meta_line unique key (first chars before the =)
        unsigned int meta_line_index = meta_line.find("=", 0);
        string meta_line_prefix = meta_line.substr(0, meta_line_index);

        // check if the meta_line_prefix is in the header_lines, if so add it to the appropirate list
        if (this->header_lines.find(meta_line_prefix) != header_lines.end()) // the meta_line is a header line so replace what was there
        {
            this->header_lines[meta_line_prefix] = meta_line;
        }
        else if (header_lists.find(meta_line_prefix) != header_lists.end() &&
            !metaInfoIdExistsInVector(meta_line, this->header_lists[meta_line_prefix])) // check if the metalineprefix is in the headerLists, if so add it to the appropirate list
        {
            this->header_lists[meta_line_prefix].push_back(meta_line);
        }
    }

    string VCFHeader::getHeaderString()
    {
        // getHeaderString generates the string each time it is called
        string header_string;

        // start by adding the header_lines
        for (vector<string>::const_iterator header_lines_iter = this->header_line_names_ordered.begin(); header_lines_iter != this->header_line_names_ordered.end(); ++header_lines_iter)
        {
            if (this->header_lines[(*header_lines_iter)] != "")
            {
                header_string += this->header_lines[(*header_lines_iter)] + "\n";
            }
        }

        // next add header_lists
        for (vector<string>::const_iterator header_lists_iter = this->header_list_names_ordered.begin(); header_lists_iter != this->header_list_names_ordered.end(); ++header_lists_iter)
        {
            vector<string> tmp_header_lists = this->header_lists[(*header_lists_iter)];
            for (vector<string>::const_iterator header_list = tmp_header_lists.begin(); header_list != tmp_header_lists.end(); ++header_list)
            {
                header_string += (*header_list) + "\n";
            }
        }

        // last add header columns
        vector<string>::const_iterator last_element = this->header_columns.end() - 1;
        for (vector<string>::const_iterator header_column_iter = this->header_columns.begin(); header_column_iter != this->header_columns.end(); ++header_column_iter)
        {
            string delimiter = (header_column_iter == last_element) ? "\n" : "\t";
            header_string += (*header_column_iter) + delimiter;
        }
        return header_string;
    }

    bool VCFHeader::metaInfoIdExistsInVector(const string& meta_line, vector<string>& meta_lines)
    {
        // extract the id from meta_line
        size_t meta_line_id_start_idx = meta_line.find("ID=", 0); // used for the start of the substring index
        size_t meta_line_id_end_idx = meta_line.find(",", meta_line_id_start_idx); // used for end of the substring index
        string meta_line_id = (meta_line_id_start_idx < meta_line_id_end_idx) ? meta_line.substr(meta_line_id_start_idx, meta_line_id_end_idx - meta_line_id_start_idx) : "";

        for (vector<string>::const_iterator iter = meta_lines.begin(); iter != meta_lines.end(); ++iter)
        {
            // extract the id from iter's meta_line string
            size_t iter_meta_line_id_start_idx = (*iter).find("ID=", 0);
            size_t iter_meta_line_id_end_idx = (*iter).find(",", iter_meta_line_id_start_idx);
            string iter_meta_line_id = (iter_meta_line_id_start_idx < iter_meta_line_id_end_idx) ? (*iter).substr(iter_meta_line_id_start_idx, iter_meta_line_id_end_idx - iter_meta_line_id_start_idx) : "";
            // compare the meta_line_id with the iter_meta_line_id
            if (strcasecmp(meta_line_id.c_str(), iter_meta_line_id.c_str()) == 0)
            {
                return true;
            }
        }
        return false;
    }

    void VCFHeader::addHeaderColumn(const string& header_column)
    {
        // don't add duplicates
        //  vector<string>::iterator test = find(this->header_columns.begin(), this->header_columns.end(), header_column);
        if (find(this->header_columns.begin(), this->header_columns.end(), header_column) == this->header_columns.end())
        {
            this->header_columns.push_back(header_column);
        }
    }

} // end namespace vcf
