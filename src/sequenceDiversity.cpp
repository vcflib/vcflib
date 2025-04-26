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
#include "gpatInfo.hpp"
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


void printHelp(void){
  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << "      The sequenceDiversity program calculates two popular metrics of  haplotype diversity: pi and                                  " << endl;
  cerr << "      extended haplotype homozygoisty (eHH).  Pi is calculated using the Nei and Li 1979 formulation.                               " << endl;
  cerr << "      eHH a convenient way to think about haplotype diversity.  When eHH = 0 all haplotypes in the window                           " << endl;
  cerr << "      are unique and when eHH = 1 all haplotypes in the window are identical.                           " << endl;

  cerr << endl;
  cerr << "Output : 5 columns:"           << endl;
  cerr << "         1.  seqid"            << endl;
  cerr << "         2.  start of window"  << endl;
  cerr << "         3.  end of window  "  << endl;
  cerr << "         4.  pi             "  << endl;
  cerr << "         5.  eHH            "  << endl;
  cerr << endl << endl;
  cerr << "INFO: usage: sequenceDiversity --target 0,1,2,3,4,5,6,7 --file my.vcf                                                                      " << endl;
  cerr << endl;
  cerr << "INFO: required: t,target     -- argument: a zero base comma separated list of target individuals corresponding to VCF columns        " << endl;
  cerr << "INFO: required: f,file       -- argument: a properly formatted phased VCF file                                                       " << endl;
  cerr << "INFO: required: y,type       -- argument: type of genotype likelihood: PL, GL or GP                                                  " << endl;
  cerr << "INFO: optional: a,af         -- sites less than af  are filtered out; default is 0                                          " << endl;
  cerr << "INFO: optional: r,region     -- argument: a tabix compliant region : \"seqid:0-100\" or \"seqid\"                                    " << endl;
  cerr << "INFO: optional: w,window     -- argument: the number of SNPs per window; default is 20                                               " << endl;
  cerr << endl << "Type: statistics" << endl << endl;
  cerr << endl;

  printVersion();

  exit(1);
}

void clearHaplotypes(std::vector<std::pair<std::string, std::string>>& haplotypes){
  for(int i= 0; i < haplotypes.size(); i++){
    haplotypes[i].first.clear();
    haplotypes[i].second.clear();
  }
}

void pi(map<string, int> & hapWin, int nHaps, double * pi, double * eHH, int wlen){

  double nchooseSum = 0;
  // summing over all possible haplotypes
  for(const auto& it : hapWin){
    nchooseSum += r8_choose(it.second, 2);
  }

  double piSum = 0;
  // all unique pairwise
  for(map<string, int>::iterator it = hapWin.begin();
      it != hapWin.end(); it++){

    // advancing it
    map<string, int>::iterator iz = it;
    iz++;
    for(;iz != hapWin.end(); iz++){
      // different bases
      int ndiff = 0;
      for(int i = 0; i < it->first.size();i++){
	if(it->first[i] != iz->first[i]){
	  ndiff += 1;
	}
      }
      double f1 = double(it->second)/double(nHaps);
      double f2 = double(iz->second)/double(nHaps);
      double perBaseDiff = double(ndiff)/double(wlen);

      piSum += f1*f2*perBaseDiff;
    }
  }


  *pi  = piSum;
  *eHH = nchooseSum / r8_choose(nHaps, 2);
}


//calc(haplotypes, nsamples, positions, targetAFS, backgroundAFS, external, derived, windowSize, target_h, background_h, currentSeqid)
void calc(const std::vector<std::pair<std::string, std::string>>& haplotypes, int, 
    vector<long int> pos, const vector<double>&, const vector<double>&, int ,
    int, int window,  const vector<int> & target, const vector<int> & , const string& seqid){

  if(haplotypes[0].first.length() < (window-1) ){
    return;
  }

  for(int snpA = 0; snpA < haplotypes[0].first.length() - window; snpA += 1){

    map <string, int> targetHaplotypes;


    for(const auto t : target){

      string haplotypeA;
      string haplotypeB;

      haplotypeA += haplotypes[t].first.substr(snpA, window) ;
      haplotypeB += haplotypes[t].second.substr(snpA, window) ;

      targetHaplotypes[haplotypeA]++;
      targetHaplotypes[haplotypeB]++;

    }

    double piEst;
    double eHH = 0;

    // target haplotype are the counts of the unique haplotypes

    int wlen = pos[snpA + window] - pos[snpA];

    pi(targetHaplotypes, target.size()*2, &piEst, &eHH, wlen);

    cout << seqid << "\t" << pos[snpA] << "\t" << pos[snpA + window] << "\t" << piEst << "\t" << eHH << endl;

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

  map<int, int> targetIndex, backgroundIndex;

  // deltaaf is the difference of allele frequency we bother to look at

  // ancestral state is set to zero by default


  int counts = 0;

  // phased

  int phased   = 0;

  // use the background allele frequency

  int external = 0;

  // "11" vs "00"

  int derived = 0;

  int windowSize = 20;

  // allele frequency to filter out
  double af_filt = 0;

  string type = "NA";

    const struct option longopts[] =
      {
	{"version"     , 0, 0, 'v'},
	{"help"        , 0, 0, 'h'},
        {"file"        , 1, 0, 'f'},
	{"target"      , 1, 0, 't'},
	{"region"      , 1, 0, 'r'},
	{"type"        , 1, 0, 'y'},
	{"window"      , 1, 0, 'w'},
	{"external"    , 1, 0, 'e'},
	{"af"          , 1, 0, 'a'},
	{"derived"     , 1, 0, 'd'},
	{0,0,0,0}
      };

    int findex;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "a:w:y:r:t:b:f:edhv", longopts, &findex);

	switch (iarg)
	  {
	  case 'h':
	    {
	      printHelp();
	      break;
	    }
	  case 'v':
	    {
	      printVersion();
	      break;
	    }
	  case 'y':
	    {
	      type = optarg;
	      break;
	    }
	  case 'a':
	    {
	      af_filt = atof(optarg);
	      cerr << "INFO: filtering out allele frequencies less than: " << af_filt << endl;
	      break;
	    }
	  case 't':
	    {
	      loadIndices(targetIndex, optarg);
	      cerr << "INFO: there are " << targetIndex.size() << " individuals in the target" << endl;
	      cerr << "INFO: target ids: " << optarg << endl;
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
	  case 'e':
	    {
	      external = 1;
	      cerr << "INFO: using background group\'s allele frequecy" << endl;
	      break;
	    }
	  case 'd':
	    {
	      derived == 1;
	      cerr << "INFO: count haplotypes \"11\" rather than \"00\"" << endl;
	      break;
	    }
	  case 'w':
	    {
	      string win = optarg;
	      windowSize = atof( win.c_str() );
	      break;
	    }
	  default :
	    break;
	  }
      }

    map<string, int> okayGenotypeLikelihoods;
    okayGenotypeLikelihoods["PL"] = 1;
    okayGenotypeLikelihoods["GL"] = 1;
    okayGenotypeLikelihoods["GP"] = 1;
    okayGenotypeLikelihoods["GT"] = 1;

    if(targetIndex.size() < 2){
      cerr << endl;
      cerr << "FATAL: failed to specify a target - or - too few individuals in the target" << endl;
      printHelp();
      return 1;
    }

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

    cerr << "INFO: window size: " << windowSize << endl;

    variantFile.open(filename);

   if(region != "NA"){
     variantFile.setRegion(region);
   }

    if (!variantFile.is_open()) {
        return 1;
    }

    Variant var(variantFile);

    vector<string> samples = variantFile.sampleNames;
    int nsamples = samples.size();

    vector<int> target_h, background_h;

    int index = 0, indexi = 0;

    for(const auto& _ : samples){
      if(targetIndex.find(index) != targetIndex.end() ){
	    target_h.push_back(indexi);
	    indexi++;
      }
      if(backgroundIndex.find(index) != backgroundIndex.end()){
	background_h.push_back(indexi);
	indexi++;
      }
      index++;
    }

    vector<long int> positions;
    vector<double>   targetAFS;
    vector<double>   backgroundAFS;

    std::vector<std::pair<std::string, std::string>> haplotypes(nsamples);
 
    string currentSeqid = "NA";

    while (variantFile.getNextVariant(var)) {

      if(!var.isPhased()){
	cerr <<"FATAL: Found an unphased variant. All genotypes must be phased!" << endl;
	printHelp();
	return(1);
      }
      if(var.alt.size() > 1){
	continue;
      }
      if(currentSeqid != var.sequenceName){
	if(haplotypes[0].first.length() > windowSize){
	  calc(haplotypes, nsamples, positions, targetAFS, backgroundAFS, external, derived, windowSize, target_h, background_h, currentSeqid);
	}
	clearHaplotypes(haplotypes);
	positions.clear();
	currentSeqid = var.sequenceName;
	targetAFS.clear();
	backgroundAFS.clear();
      }


      vector < map< string, vector<string> > > target, background, total;

      int sindex = 0;

      for(int nsamp = 0; nsamp < nsamples; nsamp++){

	map<string, vector<string> > sample = var.samples[ samples[nsamp]];

	if(targetIndex.find(sindex) != targetIndex.end() ){
	  target.push_back(sample);
	  total.push_back(sample);
	}
	if(backgroundIndex.find(sindex) != backgroundIndex.end()){
	  background.push_back(sample);
	  total.push_back(sample);
	}
	sindex += 1;
      }

      std::unique_ptr<genotype> populationTarget    ;
      std::unique_ptr<genotype> populationBackground;
      std::unique_ptr<genotype> populationTotal     ;


      if(type == "PL"){
	populationTarget     = std::make_unique<pl>();
	populationBackground = std::make_unique<pl>();
	populationTotal      = std::make_unique<pl>();
      }
      if(type == "GL"){
	populationTarget     = std::make_unique<gl>();
	populationBackground = std::make_unique<gl>();
	populationTotal      = std::make_unique<gl>();
      }
      if(type == "GP"){
	populationTarget     = std::make_unique<gp>();
	populationBackground = std::make_unique<gp>();
	populationTotal      = std::make_unique<gp>();
      }
      if(type == "GT"){
  populationTarget     = std::make_unique<gt>();
  populationBackground = std::make_unique<gt>();
  populationTotal      = std::make_unique<gt>();
      }

      populationTarget->loadPop(target,         var.position);

      populationBackground->loadPop(background, var.position);

      populationTotal->loadPop(total,           var.position);

      if(populationTotal->af < af_filt){
	continue;
      }

      targetAFS.push_back(populationTarget->af);
      backgroundAFS.push_back(populationBackground->af);
      positions.push_back(var.position);
      loadPhased(haplotypes, populationTotal.get());

    }

    calc(haplotypes, nsamples, positions, targetAFS, backgroundAFS, external, derived, windowSize, target_h, background_h, currentSeqid);

    return 0;
}
