/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"
#include <getopt.h>

using namespace std;
using namespace vcflib;

void printSummary(char** argv) {
  cerr << "usage: vcf2tsv [-n null_string] [-g]" << " [vcf file]" << endl << endl
       << "Converts VCF to per-allelle or per-genotype tab-delimited format, using null string to replace empty values in the table." << endl
       << "Specifying -g will output one line per sample with genotype information." << endl
       << "When there is more than one alt allele there will be multiple rows, one for each allele and, the info will match the 'A' index" << endl;
    cerr << endl << "Type: transformation" << endl << endl;
    exit(1);
}

void loadGenoSS(std::stringstream & ss,
                std::stringstream & info,
                vcflib::Variant & var,
                vcflib::VariantCallFile & vcf,
                const std::string & nullval,
                const std::vector<std::string> & formatfields){
    for (const auto& s : var.samples) {
        const string& sampleName = s.first;
        ss << info.str() << "\t" << sampleName;
        const map<string, vector<string> >& sample = s.second;
        for (const auto& f : formatfields) {
            const auto sampleIt = sample.find(f);
        	if (sampleIt != sample.end()) {
                ss << "\t" << join(sampleIt->second, ",");
            } else {
                ss << "\t" << nullval;
            }
        }
        ss << std::endl;
    }
}

void loadInfoSS(std::stringstream & ss,
                std::map<std::string, bool> & infoToKeep,
                vcflib::Variant & var,
                vcflib::VariantCallFile & vcf,
                std::string & nullval,
                std::vector<std::string> & formatfields,
                bool genotype){

    int altindex = 0;

    for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a, ++altindex) {
        string& altallele = *a;

        std::stringstream o;
        o   << var.sequenceName << "\t"
            << var.position << "\t"
            << var.id << "\t"
            << var.ref << "\t"
            << altallele << "\t"
            << var.quality << "\t"
            << var.filter;


        for(const auto& it : infoToKeep) {
        	if(var.info.find(it.first) == var.info.end()){
                o << "\t" << nullval ;
            }
            else{
                if(it.second == true){
                    o << "\t" << var.info[it.first].front();
                }
                else{
                    if(vcf.infoCounts.find(it.first) == vcf.infoCounts.end()){

                        if(vcf.infoCounts[it.first] == ALLELE_NUMBER){
                            o << "\t" << var.info[it.first].at(altindex) ;
                        }
                    }
                    else{
                        o << "\t" << join(var.info[it.first], ",") ;
                    }
                }
            }
        }
        if(genotype){
            loadGenoSS(ss, o, var, vcf, nullval, formatfields);
        }
        else{
            ss << o.str() << std::endl;
        }
    }
}

int main(int argc, char** argv) {

    string nullval = ".";
    bool genotypes = false;

    int c;
    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
            {"null-value", required_argument, 0, 'n'},
            {"genotypes", no_argument, 0, 'g'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hn:g",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {

	    case 'n':
	        nullval = optarg;
            break;

        case 'g':
            genotypes = true;
            break;

        case 'h':
            printSummary(argv);
            break;

        case '?':
            printSummary(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    VariantCallFile variantFile;
    bool usingstdin = false;
    string inputFilename;
    if (optind == argc - 1) {
        inputFilename = argv[optind];
        variantFile.open(inputFilename);
    } else {
        if (!variantFile.open(std::cin)) {
            if (argc == 1) {
                printSummary(argv);
            } else {
                cerr << "could not open stdin for reading as VCF" << endl;
                exit(1);
            }
        }
        usingstdin = true;
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    // obtain all possible field names
    // true means it a bool field flag
    std::map<std::string, bool> keepFields;

    for (map<string, VariantFieldType>::iterator i = variantFile.infoTypes.begin(); i != variantFile.infoTypes.end(); ++i) {
        if (i->second == FIELD_BOOL) {
            keepFields[i->first] = true;
        } else {
            keepFields[i->first] = false;
        }
    }
    vector<string> formatfields;
    if (genotypes) {
        for (map<string, VariantFieldType>::iterator f = variantFile.formatTypes.begin(); f != variantFile.formatTypes.end(); ++f) {
            formatfields.push_back(f->first);
        }
    }

    // write header
    // defaults
    std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER";

    for (std::map<std::string, bool>::iterator i = keepFields.begin(); i != keepFields.end(); ++i) {
        cout << "\t" << i->first;
    }

    if (genotypes) {
        cout << "\t" << "SAMPLE";
        for (vector<string>::iterator f = formatfields.begin(); f != formatfields.end(); ++f) {
            cout << "\t" << *f;
        }
    }
    std::cout << std::endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        stringstream outputRecord;

        loadInfoSS(outputRecord, keepFields, var, variantFile, nullval, formatfields, genotypes);

        std::cout << outputRecord.str() ;

    }
    return 0;
}
