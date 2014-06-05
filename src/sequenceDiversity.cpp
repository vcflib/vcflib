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

using namespace std;
using namespace vcf;
void printVersion(void){
	    cerr << "INFO: version 1.1.0 ; date: April 2014 ; author: Zev Kronenberg; email : zev.kronenberg@utah.edu " << endl;
	    exit(1);
}

void printHelp(void){
  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << "      The sequenceDiversity program calculates two popular metrics of  haplotype diversity: pi and                                  " << endl;
  cerr << "      extended haplotype homozygoisty (eHH).  Pi is calculated using the Nei and Li 1979 formulation.                               " << endl;
  cerr << "      eHH a convenient way to think about haplotype diversity.  When eHH = 0 all haplotypes in the window                           " << endl;
  cerr << "      are unique and when eHH = 1 all haplotypes in the window are identical. The window size is 20 SNPs.                           " << endl;

  cerr << endl;
  cerr << "Output : 5 columns:"            << endl;
  cerr << "         1.  seqid"            << endl;
  cerr << "         2.  start of window"  << endl;
  cerr << "         3.  end of window  "  << endl;
  cerr << "         4.  pi             "  << endl;
  cerr << "         5.  eHH            "  << endl;
  cerr << endl << endl;
  cerr << "INFO: usage: sequenceDiversity --target 0,1,2,3,4,5,6,7 --file my.vcf                                                                      " << endl;
  cerr << endl;
  cerr << "INFO: required: t,target     -- argument: a zero based comma separated list of target individuals corrisponding to VCF columns        " << endl;
  cerr << "INFO: required: f,file       -- argument: a properly formatted phased VCF file                                                       " << endl;
  cerr << "INFO: required: y,type       -- argument: type of genotype likelihood: PL, GL or GP                                                  " << endl;
  cerr << "INFO: optional; r,region     -- argument: a tabix compliant region : \"seqid:0-100\" or \"seqid\"                                    " << endl; 
  cerr << endl;
 
  printVersion();

  exit(1);
}

void clearHaplotypes(string haplotypes[][2], int ntarget){
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

double pi(map<string, int> & hapCounts, int nHaps, double * pi, double * eHH){

  double sum = 0;

  double nchooseSum = 0;

  for(map<string, int>::iterator firstSeq = hapCounts.begin(); firstSeq != hapCounts.end(); firstSeq++){
    
    nchooseSum += r8_choose(firstSeq->second, 2);

    for(map<string, int>::iterator secondSeq = hapCounts.begin(); secondSeq != hapCounts.end(); secondSeq++){
      if(firstSeq->first == secondSeq->first){
	continue;
      }
      
      double seqDiff = 0;
      
      for(int i = 0; i < (firstSeq->first).length(); i++){
	if((firstSeq->first).substr(i, 1) == (secondSeq->first).substr(i, 1)){
	  seqDiff += 1;
	}
      }
      sum += ((double(firstSeq->second) / nHaps ) *  (double(secondSeq->second) / nHaps ) * seqDiff); 
    }
    hapCounts.erase(firstSeq->first);
  }
  
  *pi  = sum;
  *eHH = nchooseSum / r8_choose(nHaps, 2);
  
}

void calc(string haplotypes[][2], int nhaps, vector<long int> pos, vector<double> tafs, vector<double> bafs, int external, long int window, int derived, vector<int> & target, vector<int> & background, string seqid){
  
  for(int long snpA = 0; snpA < haplotypes[0][0].length() - 20; snpA += 1){
    
    map <string, int> targetHaplotypes;
    
    for(int targetIndex = 0; targetIndex < target.size(); targetIndex++ ){
      
      string haplotypeA;
      string haplotypeB;
      
      haplotypeA += haplotypes[target[targetIndex]][0].substr(snpA, 20) ;
      haplotypeB += haplotypes[target[targetIndex]][1].substr(snpA, 20) ;
      
      targetHaplotypes[haplotypeA]++;
      targetHaplotypes[haplotypeB]++;
      
    }
    
    double piEst;
    double eHH;

    pi(targetHaplotypes, target.size()*2, &piEst, &eHH); 

    cout << seqid << "\t" << pos[snpA] << "\t" << pos[snpA + 20] << "\t" << piEst << "\t" << eHH << endl;
   
  }

}

void loadPhased(string haplotypes[][2], genotype * pop, int ntarget){
  
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

    if(it.size() < 2){
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

    string haplotypes [nsamples][2];    
    
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


      map<string, map<string, vector<string> > >::iterator s     = var.samples.begin(); 
      map<string, map<string, vector<string> > >::iterator sEnd  = var.samples.end();
      
      vector < map< string, vector<string> > > target, background, total;
      
      int sindex = 0;
      
      for (; s != sEnd; s++) {	  
	
	map<string, vector<string> >& sample = s->second;
	  
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
      
      targetAFS.push_back(populationTarget->af);
      backgroundAFS.push_back(populationBackground->af);
      positions.push_back(var.position);
      loadPhased(haplotypes, populationTotal, nsamples);
      
    }

    calc(haplotypes, nsamples, positions, targetAFS, backgroundAFS, external, derived, windowSize, target_h, background_h, currentSeqid);
    
    return 0;		    
}