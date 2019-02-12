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

  cerr << "INFO: LD --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf -e -d -r                                           " << endl;
  cerr << endl;
  
  cerr << "INFO: required: t,target     -- argument: a zero base comma seperated list of target individuals corrisponding to VCF columns        " << endl;
  cerr << "INFO: required: b,background -- argument: a zero base comma seperated list of background individuals corrisponding to VCF columns    " << endl;
  cerr << "INFO: required: f,file       -- argument: a properly formatted phased VCF file                                                       " << endl;
  cerr << "INFO: required: y,type       -- argument: type of genotype likelihood: PL, GL or GP                                                  " << endl;
  cerr << "INFO: optional: w,window     -- argument: window size to average LD; default is 1000                                                 " << endl;
  cerr << "INFO: optional: e,external   -- switch: population to calculate LD expectation; default is target                                    " << endl;
  cerr << "INFO: optional: d,derived    -- switch: which haplotype to count \"00\" vs \"11\"; default \"00\",                                   " << endl;
  cerr << endl;
 
  printVersion();

  exit(1);
}

void clearHaplotypes(string **haplotypes, int ntarget){
  for(int i= 0; i < ntarget; i++){
    haplotypes[i][0].clear();
    haplotypes[i][1].clear();
  }
}

void loadIndices(map<int, int> & index, string set){
  
  vector<string>  indviduals = split(set, ",");
  vector<string>::iterator it = indviduals.begin();
  
  for(; it != indviduals.end(); it++){
    index[ atoi( (*it).c_str() ) ] = 1;
  }
}

void calc(string **haplotypes, int nhaps, vector<long int> pos, vector<double> tafs, vector<double> bafs, int external, long int window, int derived, vector<int> & target, vector<int> & background, string seqid){

  for(int long snpA = 0; snpA < haplotypes[0][0].length() - 100; snpA++){

    double sumLD = 0;
    double nLD   = 0;

    for(int snpB = snpA; snpB < (snpA + 100); snpB++){
      
      map<string, int> targetHaplotypes;
      
      targetHaplotypes["11"] = 1;
      targetHaplotypes["01"] = 1;
      targetHaplotypes["10"] = 1;
      targetHaplotypes["00"] = 1;
      
      for(int targetIndex = 0; targetIndex < target.size(); targetIndex++ ){
	
	string haplotypeA;
	string haplotypeB;
	
	haplotypeA += haplotypes[target[targetIndex]][0].substr(snpA, 1) +=  haplotypes[target[targetIndex]][0].substr(snpB, 1);
	haplotypeB += haplotypes[target[targetIndex]][1].substr(snpA, 1) +=  haplotypes[target[targetIndex]][1].substr(snpB, 1); 
	
	targetHaplotypes[haplotypeA]++;
	targetHaplotypes[haplotypeB]++;
	  
      }
        

      double sa = tafs[snpA];
      double sb = tafs[snpA];

      string hapID = "00";
      
      if(external == 1){
	sa = bafs[snpA];
	sb = bafs[snpB];
      }
      if(derived == 1){
	hapID = "11";
      }
      else{
	sa = 1 - sa;
	sb = 1 - sb;
      }

      double observation = targetHaplotypes[hapID] / ( target.size()*2 ) ;      
      double expectation = sa*sb;

      double d           = observation - expectation; 
      //double r           = d / sqrt((1 - afs[snpA])*(1-afs[snpB])*afs[snpA]*afs[snpB]);
    
      nLD   += 1;
      sumLD += d;
      
    }
    
    cout << seqid << "\t" << pos[snpA] << "\t" << tafs[snpA] << "\t" << sumLD << "\t" << nLD  << endl;

  }

}

void loadPhased(string **haplotypes, genotype * pop, int ntarget){
  
  int indIndex = 0;

  for(vector<string>::iterator ind = pop->gts.begin(); ind != pop->gts.end(); ind++){
    string g = (*ind);
    vector< string > gs = split(g, "|");
    haplotypes[indIndex][0].append(gs[0]);
    haplotypes[indIndex][1].append(gs[1]);
    indIndex += 1;
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

  int phased   = 0;
  
  // use the background allele frequency

  int external = 0;

  // "11" vs "00"

  int derived = 0;

  long int windowSize = 1000;

  string type = "NA";

    const struct option longopts[] = 
      {
	{"version"     , 0, 0, 'v'},
	{"help"        , 0, 0, 'h'},
        {"file"        , 1, 0, 'f'},
	{"target"      , 1, 0, 't'},
	{"background"  , 1, 0, 'b'},
	{"region"      , 1, 0, 'r'},
	{"type"        , 1, 0, 'y'},
	{"window"      , 1, 0, 'w'},
	{"external"    , 1, 0, 'e'},
	{"derived"     , 1, 0, 'd'},
	{0,0,0,0}
      };

    int findex;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "w:y:r:t:b:f:edhv", longopts, &findex);
	
	switch (iarg)
	  {
	  case 'h':
	    {
	      printHelp();
	    }
	  case 'v':
	    {
	      printVersion();
	    }
	  case 'y':
	    {
	      type = optarg;
	      break;
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
	      windowSize = atol( win.c_str() );
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
    
   if(region != "NA"){
     if(! variantFile.setRegion(region)){
       cerr <<"FATAL: unable to set region" << endl;
       return 1;
     }
   }
    
    if (!variantFile.is_open()) {
        return 1;
    }
    
    Variant var(variantFile);

    vector<string> samples = variantFile.sampleNames;
    int nsamples = samples.size();

    vector<int> target_h, background_h;

    int index, indexi = 0;

    for(vector<string>::iterator samp = samples.begin(); samp != samples.end(); samp++){
     
      if(it.find(index) != it.end() ){
	target_h.push_back(indexi);
	indexi++;
      }
      if(ib.find(index) != ib.end()){
	background_h.push_back(indexi);
	indexi++;
      }
      index++;
    }
    
    vector<long int> positions;
    vector<double>   targetAFS;
    vector<double>   backgroundAFS;

    string **haplotypes = new string*[nsamples];
	for (int i = 0; i < nsamples; i++) {
	  haplotypes[i] = new string[2];
	}
    
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
	if(haplotypes[0][0].length() > 10){
	  calc(haplotypes, nsamples, positions, targetAFS, backgroundAFS, external, derived, windowSize, target_h, background_h, currentSeqid);
	}
	clearHaplotypes(haplotypes, nsamples);
	positions.clear();
	currentSeqid = var.sequenceName;
	targetAFS.clear();
	backgroundAFS.clear();
      }
      
      vector < map< string, vector<string> > > target, background, total;
      
      int sindex = 0;
      
      for(int nsamp = 0; nsamp < nsamples; nsamp++){

	map<string, vector<string> > sample = var.samples[ samples[nsamp]];

	if(it.find(sindex) != it.end() ){
	  target.push_back(sample);
	  total.push_back(sample);	  
	}
	if(ib.find(sindex) != ib.end()){
	  background.push_back(sample);
	  total.push_back(sample);	  
	}	
	sindex += 1;
      }
      
      genotype * populationTarget    ;
      genotype * populationBackground;
      genotype * populationTotal     ;
      
      if(type == "PL"){
	populationTarget     = new pl();
	populationBackground = new pl();
	populationTotal      = new pl();
      }
      if(type == "GL"){
	populationTarget     = new gl();
	populationBackground = new gl();
	populationTotal      = new gl();
      }
      if(type == "GP"){
	populationTarget     = new gp();
	populationBackground = new gp();
	populationTotal      = new gp();
      }
      
      populationTarget->loadPop(target,         var.sequenceName, var.position);
      
      populationBackground->loadPop(background, var.sequenceName, var.position);
	
      populationTotal->loadPop(total,           var.sequenceName, var.position);
      
      
      if(populationTotal->af > 0.95 || populationTotal->af < 0.05){
	continue;
      }

      targetAFS.push_back(populationTarget->af);
      backgroundAFS.push_back(populationBackground->af);
      positions.push_back(var.position);
      loadPhased(haplotypes, populationTotal, nsamples);
      
    }

    calc(haplotypes, nsamples, positions, targetAFS, backgroundAFS, external, derived, windowSize, target_h, background_h, currentSeqid);
    
    return 0;		    
}
