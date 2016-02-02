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

struct indv{
  int nhet  ;
  int nhom  ;
  int nalt  ;
  int nocall;
};



void printHelp(void){

  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << "      Summarizes genotype counts for bi-allelic SNVs and indel " << endl;

  cerr << "INFO: output: table of genotype counts for each individual." << endl;


  cerr << "INFO: usage:  genotypeSummmary --type PL --target 0,1,2,3,4,5,6,7 --file my.vcf --snp                                                               " << endl;
  cerr << endl;
  cerr << "INFO: required: t,target     -- a zero based comma separated list of target individuals corresponding to VCF columns        " << endl;
  cerr << "INFO: required: f,file       -- proper formatted VCF                                                                        " << endl;
  cerr << "INFO: required, y,type       -- genotype likelihood format; genotype : GL,PL,GP                                             " << endl;
  cerr << "INFO: optional, r,region     -- a tabix compliant region : chr1:1-1000 or chr1                                              " << endl;
  cerr << "INFO: optional, s,snp        -- Only count SNPs                                              " << endl;

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

  bool snp = false;

  // set the random seed for MCMC

  srand((unsigned)time(NULL));

  // the filename

  string filename;

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
	{"snp"       , 0, 0, 's'},
	{0,0,0,0}
      };

    int index;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "y:r:d:t:b:f:chvs", longopts, &index);
	
	switch (iarg)
	  {
	  case 's':
	    {
	      snp = true;
	      break;
	    }
	  case 'h':
	    {
	      printHelp();
	      return 0;
	    }
	  case 'v':
	    {
	      printVersion();
	      return 0;
	    }
	  case 't':
	    {
	      loadIndices(it, optarg);
	      cerr << "INFO: there are " << it.size() << " individuals in the target" << endl;
	      cerr << "INFO: target ids: " << optarg << endl;
	      break;
	    }
	  case 'b':
	    {
	      loadIndices(ib, optarg);
	      cerr << "INFO: there are " << ib.size() << " individuals in the background" << endl;
	      cerr << "INFO: background ids: " << optarg << endl;
	      break;
	    }
	  case 'f':
	    {
	      cerr << "INFO: file: " << optarg  <<  endl;
	      filename = optarg;
	      break;
	    }
	  case 'r':
	    {
	      cerr << "INFO: set seqid region to : " << optarg << endl;
	      region = optarg; 
	      break;
	    }
	  case 'y':
	    {
	      type = optarg;
	      cerr << "INFO: set genotype likelihood to: " << type << endl;
	      break;
	    }
	  default:
	    break;
	  }

      }

    if(filename.empty()){
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

    vector<indv *> countData;
    

    for(int i = 0; i < it.size(); i++){
      indv * dip = new indv;

      dip->nhet   = 0;
      dip->nhom   = 0;
      dip->nalt   = 0;
      dip->nocall = 0;

      countData.push_back(dip);
    }
    
    
    while (variantFile.getNextVariant(var)) {
        
	// biallelic sites naturally 

	if(var.alt.size() > 1){
	  continue;
	}
	if(snp){
	  bool hit =false;

	  for(vector<string>::iterator it = var.alleles.begin(); it != var.alleles.end(); it++){
	    if((*it).size() > 1){
	      hit = true;
	    }
	  }
	  if(hit){
	    continue;
	  }
	}
	vector < map< string, vector<string> > > target, background, total;
	        
	int index = 0;

	for(int nsamp = 0; nsamp < nsamples; nsamp++){

	  map<string, vector<string> > sample = var.samples[ samples[nsamp]];
	      if(it.find(index) != it.end() ){
		target.push_back(sample);
	    }            
	    index += 1;
	}
	
	genotype * populationTarget      ;

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
	
	for(int i = 0; i < populationTarget->genoIndex.size() ; i++){
	  if(populationTarget->genoIndex[i] == -1){
	    countData[i]->nocall += 1;
	  }
	  else if(populationTarget->genoIndex[i] == 0){
	    countData[i]->nhom += 1;
	  }
	  else if(populationTarget->genoIndex[i] == 1){
	    countData[i]->nhet += 1;
	  }
	  else if(populationTarget->genoIndex[i] == 2){
	    countData[i]->nalt += 1;
	  }
	  else{
	    std::cerr << "FATAL: unkown genotype index" << std::endl;
	    exit(1);
	  }
	}
	delete populationTarget;

    }
    
    std::cout << "#sample-id\tn-nocall\tn-hom-ref\tn-het\tn-alt" << std::endl;
    for(int i = 0; i < countData.size(); i++){
      std::cout << samples[i] 
		<< "\t" << countData[i]->nocall 
		<< "\t" << countData[i]->nhom
		<< "\t" << countData[i]->nhet
		<< "\t" << countData[i]->nalt
		<< std::endl;
    }


    return 0;		    
}
