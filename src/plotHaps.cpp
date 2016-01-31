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
using namespace vcflib;
void printVersion(void){
	    cerr << "INFO: version 1.0.1 ; date: April 2014 ; author: Zev Kronenberg; email : zev.kronenberg@utah.edu " << endl;
	    exit(1);
}

void printHelp(void){
  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << " plotHaps provides the formatted output that can be used with \'bin/plotHaplotypes.R\'. " << endl;

  cerr << "Output : haplotype matrix and positions" << endl << endl;

  cerr << "INFO: plotHaps  --target 0,1,2,3,4,5,6,7  --file my.phased.vcf.gz                                                           " << endl << endl;
  cerr << "INFO: required: t,target     -- argument: a zero base comma separated list of target individuals corrisponding to VCF column s        " << endl;
  cerr << "INFO: required: r,region     -- argument: a tabix compliant genomic range : \"seqid:start-end\" or \"seqid\"                          " << endl;
  cerr << "INFO: required: f,file       -- argument: proper formatted phased VCF file                                                            " << endl;
  cerr << "INFO: required: y,type       -- argument: genotype likelihood format: PL,GP,GP                                                        " << endl;
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

void calc(string **haplotypes, int nhaps, vector<double> afs, vector<long int> pos, vector<int> & target, vector<int> & background, string seqid){

  for(int snp = 0; snp < haplotypes[0][0].length(); snp++){
    
    double ehhA = 1;
    double ehhR = 1;

    double iHSA = 1;
    double iHSR = 1;

    int start = snp;
    int end   = snp;
    int core  = snp; 

    while( ehhA > 0.05 && ehhR > 0.05 ) {
     
      start -= 1;
      end   += 1;
      
      if(start == -1){
	break;
      }
      if(end == haplotypes[0][0].length() - 1){
	break;
      }
      
      map<string , int> targetH;

      double sumrT = 0;
      double sumaT = 0;
      double nrefT = 0;
      double naltT = 0;

      for(int i = 0; i < nhaps; i++){
	targetH[ haplotypes[i][0].substr(start, (end - start)) ]++;
	targetH[ haplotypes[i][1].substr(start, (end - start)) ]++;
      }     
      for( map<string, int>::iterator th = targetH.begin(); th != targetH.end(); th++){    	
	if( (*th).first.substr((end-start)/2, 1) == "1"){     
	   sumaT += r8_choose(th->second, 2);  
	   naltT += th->second;
	}
	else{
	  sumrT += r8_choose(th->second, 2);  
	  nrefT += th->second;
	}
      }

      ehhA = sumaT / (r8_choose(naltT, 2));
      ehhR = sumrT / (r8_choose(nrefT, 2));

      iHSA += ehhA;
      iHSR += ehhR;
    } 
    cout << seqid << "\t" << pos[snp] << "\t" << afs[snp] << "\t" << iHSA << "\t" << iHSR << "\t" << iHSA/iHSR << endl;
  }   
}

double EHH(string **haplotypes, int nhaps){

  map<string , int> hapcounts;

  for(int i = 0; i < nhaps; i++){
    hapcounts[ haplotypes[i][0] ]++;
    hapcounts[ haplotypes[i][1] ]++;
  }

  double sum = 0;
  double nh  = 0;

  for( map<string, int>::iterator it = hapcounts.begin(); it != hapcounts.end(); it++){
    nh  += it->second; 
    sum += r8_choose(it->second, 2);
  }

  double max = (sum /  r8_choose(nh, 2));
  
  return max;

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

void printHaplotypes(string **haps, vector<int> target, vector<long int> pos){
  for(int snp = 0; snp < haps[0][1].length(); snp++){
    cout << pos[snp] << "\t" ;
    for(int ind = 0; ind < target.size(); ind++){
      cout << haps[target[ind]][0].substr(snp , 1) << "\t";
      cout << haps[target[ind]][1].substr(snp , 1) << "\t";
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
      cerr << "FATAL: must specifiy a region" << endl;
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

    string **haplotypes = new string*[target_h.size()];
	for (int i = 0; i < target_h.size(); i++) {
	  haplotypes[i] = new string[2];
	}
    
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
      
      genotype * populationTarget    ;

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
      
      positions.push_back(var.position);
      afs.push_back(populationTarget->af);
      loadPhased(haplotypes, populationTarget, populationTarget->gts.size()); 
    }
    
    printHaplotypes( haplotypes, target_h, positions);
    
    return 0;		    
}
