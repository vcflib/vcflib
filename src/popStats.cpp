#include "Variant.h"
#include "split.h"
#include "cdflib.hpp"
#include "pdflib.hpp"
#include "var.hpp"

#include <string>
#include <iostream>
#include <math.h>  
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <getopt.h>
#include "gpatInfo.hpp"

using namespace std;
using namespace vcflib;


void printHelp(void){

  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << "      General population genetic statistics for each SNP                                                                    " << endl << endl;

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

void loadIndices(map<int, int> & index, string set){
  
  vector<string>  indviduals = split(set, ",");

  vector<string>::iterator it = indviduals.begin();
  
  for(; it != indviduals.end(); it++){
    index[ atoi( (*it).c_str() ) ] = 1;
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
	
	genotype * populationTarget      ;
	genotype * populationBackground  ;

	if(type == "PL"){
	  populationTarget     = new pl();
	}
	if(type == "GL"){
	  populationTarget     = new gl();
	}
	if(type == "GP"){
	  populationTarget     = new gp();
	}
	if(type == "GT"){
          populationTarget     = new gt();
	}
	
	populationTarget->loadPop(target, var.sequenceName, var.position);

	 //cerr << "     3. target allele frequency      "    << endl;
	 //cerr << "     4. expected heterozygosity      "    << endl;
	 //cerr << "     5. observed heterozygosity      "    << endl;
	 //cerr << "     6. number of hets               "    << endl;
	 //cerr << "     7. number of homozygous ref     "    << endl;
	 //cerr << "     8. number of homozygous alt     "    << endl;
	 //cerr << "     9. target Fis                   "    << endl;

	if(populationTarget->af == -1){
	  delete populationTarget;
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

	delete populationTarget;

    }
    return 0;		    
}
