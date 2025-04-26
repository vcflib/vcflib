/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"
#include "split.h"
#include "cdflib.hpp"
#include "pdflib.hpp"
#include "var.hpp"
#include "index.hpp"
#include "phase.hpp"

#include <string>
#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <getopt.h>
#include <memory>

using namespace std;
using namespace vcflib;
void printVersion(void){
	    cerr << "INFO: version 1.0.1 ; date: April 2014 ; author: Zev Kronenberg; email : zev.kronenberg@utah.edu " << endl;
	    exit(1);
}

void printHelp(void){
  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << " plotHaps provides the formatted output that can be used with \'bin/plotHaplotypes.R\'. " << endl << endl ;

  cerr << "Output : haplotype matrix and positions" << endl << endl;

  cerr << "INFO: plotHaps  --target 0,1,2,3,4,5,6,7  --file my.phased.vcf.gz                                                           " << endl << endl;
  cerr << "INFO: required: t,target     -- argument: a zero base comma separated list of target individuals corrisponding to VCF column s        " << endl;
  cerr << "INFO: required: r,region     -- argument: a tabix compliant genomic range : \"seqid:start-end\" or \"seqid\"                          " << endl;
  cerr << "INFO: required: f,file       -- argument: proper formatted phased VCF file                                                            " << endl;
  cerr << "INFO: required: y,type       -- argument: genotype likelihood format: PL,GP,GP                                                        " << endl;
  cerr << endl << "Type: statistics" << endl << endl;
  cerr << endl;

  exit(1);
}

void printHaplotypes(const std::vector<std::pair<std::string, std::string>>& haps, const std::vector<int>& target, vector<long int> pos){
  for(int snp = 0; snp < haps[0].second.length(); snp++){
    cout << pos[snp] << "\t" ;
    for (const int t : target)
    {
      cout << haps[t].first.substr(snp , 1) << "\t";
      cout << haps[t].second.substr(snp , 1) << "\t";
    }
    cout << endl;
  }
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

  // deltaaf is the difference of allele frequency we bother to look at

  // ancestral state is set to zero by default


  int counts = 0;

  // phased

  int phased = 0;

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

    int findex;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "y:r:t:f:hv", longopts, &findex);

	switch (iarg)
	  {
	  case 'h':
	    printHelp();
	  case 'v':
	    printVersion();
	  case 'y':
	    type = optarg;
	    break;
	  case 't':
	    loadIndices(it, optarg);
	    cerr << "INFO: there are " << it.size() << " individuals in the target" << endl;
	    cerr << "INFO: target ids: " << optarg << endl;
	    break;
	  case 'f':
	    cerr << "INFO: file: " << optarg  <<  endl;
	    filename = optarg;
	    break;
	  case 'r':
            cerr << "INFO: set seqid region to : " << optarg << endl;
	    region = optarg;
	    break;
	  default:
	    break;
	  }
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

    if(filename == "NA"){
      cerr << "FATAL: did not specify a file" << endl;
      printHelp();
      return(1);
    }

    variantFile.open(filename);

    if (!variantFile.is_open()) {
      cerr << "FATAL: could not open VCF for reading" << endl;
      printHelp();
      return 1;
    }

    if(region != "NA"){
      if(! variantFile.setRegion(region)){
	cerr <<"FATAL: unable to set region" << endl;
	return 1;
      }
    }
    else{
      cerr << "FATAL: must specify a region" << endl;
      printHelp();
      return 1;
    }

    Variant var(variantFile);

    vector<string> samples = variantFile.sampleNames;
    int nsamples = samples.size();

    vector<int> target_h, background_h;

    int index  = 0;
    int indexi = 0;

    for(vector<string>::iterator samp = samples.begin(); samp != samples.end(); samp++){

      if(it.find(index) != it.end() ){
	target_h.push_back(indexi);
	indexi++;
      }
      index++;
    }

    vector<long int> positions;

    vector<double> afs;

    std::vector<std::pair<std::string, std::string>> haplotypes(target_h.size());

    string currentSeqid = "NA";


    while (variantFile.getNextVariant(var)) {

      if(!var.isPhased()){
	cerr << "FATAL: Found an unphased variant. All genotypes must be phased!" << endl;
	return(1);
      }

      if(var.alt.size() > 1){
	continue;
      }


      vector < map< string, vector<string> > > target, background, total;

      int sindex = 0;

      for(int nsamp = 0; nsamp < nsamples; nsamp++){

	map<string, vector<string> > sample = var.samples[ samples[nsamp]];

	if(it.find(sindex) != it.end() ){
	  target.push_back(sample);
	}
	sindex += 1;
      }

      std::unique_ptr<genotype> populationTarget    ;

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

      positions.push_back(var.position);
      afs.push_back(populationTarget->af);
      loadPhased(haplotypes, populationTarget.get());
    }

    printHaplotypes( haplotypes, target_h, positions);

    return 0;
}
