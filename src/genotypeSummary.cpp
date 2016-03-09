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
  cerr << "INFO: optional, a,ancestral  -- describe counts relative to the ancestral allele defined as AA in INFO" << endl;

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

  // open standardout 

  // set region to scaffold

  string region = "NA";

  // using vcflib; thanks to Erik Garrison

 VariantCallFile variantFile;

  // zero based index for the target and background indivudals

  map<int, int> it, ib;

  // genotype likelihood format

  string type = "NA";

  // are we polarizing the counts relative to the ancestral allele?
  bool use_ancestral_state = false;
  set<char> allowed_ancestral_bases = { 'A', 'T', 'C', 'G' };

    const struct option longopts[] =
      {
	{"version"   , 0, 0, 'v'},
	{"help"      , 0, 0, 'h'},
        {"file"      , 1, 0, 'f'},
	{"target"    , 1, 0, 't'},
	{"region"    , 1, 0, 'r'},
	{"type"      , 1, 0, 'y'},
	{"snp"       , 0, 0, 's'},
	{"ancestral" , 0, 0, 'a'},
	{0,0,0,0}
      };

    int index;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "y:r:d:t:b:f:chvsa", longopts, &index);

	switch (iarg)
	  {
    case 'a':
      {
        use_ancestral_state = true;
        break;
      }
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

    bool is_open; 

    if (filename == "-") {

        is_open=variantFile.open(std::cin);

    } else {

    	is_open=variantFile.open(filename); 
	
     }
    
    if (!is_open)  {
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
    vector<string > countDataSampleName;

    for ( map<int ,int>::iterator x=it.begin(); x!=it.end(); ++x) {

        countDataSampleName.push_back(samples[x->first] ); 
    }


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

  // decide if we can polarize the site if we are using the ancestral allele
  bool ref_is_ancestral_allele = true;
  if (use_ancestral_state) {
    // we need the ancestral allele to decide what to do at this site
    if (var.info.find("AA") == var.info.end()) continue;
    string ancestral_allele = var.info["AA"].front();
    // if we do not have a polarized site with only allowed bases in the ancestral allele, skip it
    bool allowed = true;
    for (string::iterator c = ancestral_allele.begin(); c != ancestral_allele.end(); ++c) {
      if (!allowed_ancestral_bases.count(*c)) {
        allowed = false;
        break;
      }
    }
    if (!allowed) continue;
    ref_is_ancestral_allele = (ancestral_allele == var.ref);
  }

	vector < map< string, vector<string> > > target, background, total;

	int index = 0;

	for(int nsamp = 0; nsamp < nsamples; nsamp++){

	    if(it.find(index) != it.end() ){
	        const map<string, vector<string> >& sample = var.samples[ samples[nsamp]];
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
	  else if (populationTarget->genoIndex[i] == 0) {
            if (!use_ancestral_state || ref_is_ancestral_allele) {
	      countData[i]->nhom += 1;
            } else {
	      countData[i]->nalt += 1;
            }
	  }
	  else if (populationTarget->genoIndex[i] == 1){
	    countData[i]->nhet += 1;
	  }
	  else if (populationTarget->genoIndex[i] == 2) {
            if (!use_ancestral_state || ref_is_ancestral_allele) {
	      countData[i]->nalt += 1;
            } else {
	      countData[i]->nhom += 1;
            }
	  }
	  else{
	    std::cerr << "FATAL: unkown genotype index" << std::endl;
cerr << populationTarget->genoIndex[i] << endl;
cerr << var << endl;
	    exit(1);
	  }
	}
	delete populationTarget;

    }

    if (!use_ancestral_state) {
        std::cout << "#sample-id\tn-nocall\tn-hom-ref\tn-het\tn-hom-alt" << std::endl;
    } else {
        std::cout << "#sample-id\tn-nocall\tn-hom-ancestral\tn-het\tn-hom-derived" << std::endl;
    }
    for(int i = 0; i < countData.size(); i++){
        std::cout << countDataSampleName[i]
                  << "\t" << countData[i]->nocall
                  << "\t" << countData[i]->nhom
                  << "\t" << countData[i]->nhet
                  << "\t" << countData[i]->nalt
                  << std::endl;
    }


    return 0;
}
