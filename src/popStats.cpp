/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"
#include "var.hpp"
#include "index.hpp"

#include <string>
#include <iostream>
#include <cmath>
#include <ctime>
#include <getopt.h>
#include <memory>

#include "gpatInfo.hpp"

using namespace std;
using namespace vcflib;


void printHelp(void){

  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << "      General population genetic statistics for each SNP                                                                    " << endl << endl;

  cerr << R"(
  Calculates basic population statistics at bi-allelic sites. The allele frequency is the number of non-reference alleles divided by the total number of alleles.  The expected hetrozygosity is 2*p*q, where p is the non-reference allele frequency and q is 1-p.  The observed heterozgosity is the fraction of 0/1 genotypes out of all genotypes.  The inbreeding coefficient, Fis, is the relative heterozygosity of each individual vs. compared to the target group. )" << endl << endl;

  cerr << "Output : 9 columns :                 "    << endl;
  cerr << "     1. seqid                        "    << endl;
  cerr << "     2. position                     "    << endl;
  cerr << "     3. target allele frequency      "    << endl;
  cerr << "     4. expected heterozygosity      "    << endl;
  cerr << "     5. observed heterozygosity      "    << endl;
  cerr << "     6. number of hets               "    << endl;
  cerr << "     7. number of homozygous ref     "    << endl;
  cerr << "     8. number of homozygous alt     "    << endl;
  cerr << "     9. target Fis                   "    << endl;


  cerr << "INFO: usage:  popStat --type PL --target 0,1,2,3,4,5,6,7 --file my.vcf                                                                " << endl;
  cerr << endl;
  cerr << "INFO: required: t,target     -- a zero based comma separated list of target individuals corresponding to VCF columns        " << endl;
  cerr << "INFO: required: f,file       -- proper formatted VCF                                                                        " << endl;
  cerr << "INFO: required, y,type       -- genotype likelihood format; genotype : GL,PL,GP                                             " << endl;
  cerr << "INFO: optional, r,region     -- a tabix compliant region : chr1:1-1000 or chr1                                              " << endl;
  cerr << endl << "Type: statistics" << endl << endl;
  printVersion();
}


double bound(double v){
  if(v <= 0.00001){
    return  0.00001;
  }
  if(v >= 0.99999){
    return 0.99999;
  }
  return v;
}

int main(int argc, char** argv) {

  // set the random seed for MCMC

  srand((unsigned)time(NULL));

  // the filename

  string filename = "NA";

  // set region to scaffold

  string region = "NA";

  // using vcflib; thanks to Erik Garrison

  VariantCallFile variantFile;

  // zero based index for the target and background indivudals

  map<int, int> it, ib;

  // genotype likelihood format

  string type = "NA";

    const struct option longopts[] =
      {
	{"version"   , 0, 0, 'v'},
	{"help"      , 0, 0, 'h'},
        {"file"      , 1, 0, 'f'},
	{"target"    , 1, 0, 't'},
	{"region"    , 1, 0, 'r'},
	{"type"      , 1, 0, 'y'},
	{0,0,0,0}
      };

    int index;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "y:r:d:t:b:f:chv", longopts, &index);

	switch (iarg)
	  {
	  case 'h':
	    printHelp();
	    return 0;
	  case 'v':
	    printVersion();
	    return 0;
	  case 't':
	    loadIndices(it, optarg);
	    cerr << "INFO: there are " << it.size() << " individuals in the target" << endl;
	    cerr << "INFO: target ids: " << optarg << endl;
	    break;
	  case 'b':
	    loadIndices(ib, optarg);
	    cerr << "INFO: there are " << ib.size() << " individuals in the background" << endl;
	    cerr << "INFO: background ids: " << optarg << endl;
	    break;
	  case 'f':
	    cerr << "INFO: file: " << optarg  <<  endl;
	    filename = optarg;
	    break;
	  case 'r':
            cerr << "INFO: set seqid region to : " << optarg << endl;
	    region = optarg;
	    break;
	  case 'y':
	    type = optarg;
	    cerr << "INFO: set genotype likelihood to: " << type << endl;
	    break;
	  default:
	    break;
	  }

      }

    if(filename == "NA"){
      cerr << "FATAL: failed to specify a file" << endl;
      printHelp();
    }

    if(!variantFile.open(filename)){
      cerr << "FATAL: could not open file for reading" << endl;
      printHelp();
    }

    if(region != "NA"){
      if(! variantFile.setRegion(region)){
	cerr <<"FATAL: unable to set region" << endl;
	return 1;
      }
    }

    if (!variantFile.is_open()) {
      cerr << "FATAL: could not open VCF for reading" << endl;
      printHelp();
      return 1;
    }

    map<string, int> okayGenotypeLikelihoods;
    okayGenotypeLikelihoods["PL"] = 1;
    okayGenotypeLikelihoods["GL"] = 1;
    okayGenotypeLikelihoods["GP"] = 1;
    okayGenotypeLikelihoods["GT"] = 1;

    if(type == "NA"){
      cerr << "FATAL: failed to specify genotype likelihood format : PL or GL" << endl;
      printHelp();
      return 1;
    }
    if(okayGenotypeLikelihoods.find(type) == okayGenotypeLikelihoods.end()){
      cerr << "FATAL: genotype likelihood is incorrectly formatted, only use: PL or GL" << endl;
      printHelp();
      return 1;
    }

    Variant var(variantFile);

    vector<string> samples = variantFile.sampleNames;
    int nsamples = samples.size();

    while (variantFile.getNextVariant(var)) {

	// biallelic sites naturally

	if(var.alt.size() > 1){
	  continue;
	}

	vector < map< string, vector<string> > > target, background, total;

	int index = 0;

	for(int nsamp = 0; nsamp < nsamples; nsamp++){

	  map<string, vector<string> > sample = var.samples[ samples[nsamp]];

	    if(sample["GT"].front() != "./."){
	      if(it.find(index) != it.end() ){
		target.push_back(sample);
	      }
	    }
	    index += 1;
	}


	std::unique_ptr<genotype> populationTarget      ;
	unique_ptr<genotype> populationBackground  ;

	if(type == "PL"){
	  populationTarget     = std::make_unique<pl>();
	}
	if(type == "GL"){
	  populationTarget     = std::make_unique<gl>();
	}
	if(type == "GP"){
	  populationTarget     = std::make_unique<gp>();
	}
	if(type == "GT"){
    populationTarget     = std::make_unique<gt>();
	}

	populationTarget->loadPop(target, var.position);

	 //cerr << "     3. target allele frequency      "    << endl;
	 //cerr << "     4. expected heterozygosity      "    << endl;
	 //cerr << "     5. observed heterozygosity      "    << endl;
	 //cerr << "     6. number of hets               "    << endl;
	 //cerr << "     7. number of homozygous ref     "    << endl;
	 //cerr << "     8. number of homozygous alt     "    << endl;
	 //cerr << "     9. target Fis                   "    << endl;

	if(populationTarget->af == -1){
	  continue;
	}

	double ehet = 2*(populationTarget->af * (1 - populationTarget->af));

	cout << var.sequenceName << "\t"  << var.position << "\t"
	     << populationTarget->af  << "\t"
	     << ehet << "\t"
	     << populationTarget->hfrq  << "\t"
	     << populationTarget->nhet  << "\t"
	     << populationTarget->nhomr << "\t"
	     << populationTarget->nhoma << "\t"
	     << populationTarget->fis   << endl;
    }
    return 0;
}
