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
  cerr << "     pFst is a probabilistic approach for detecting differences in allele frequencies between two populations." << endl << endl;
 

  cerr << "Output : 3 columns :     "    << endl;
  cerr << "     1. seqid            "    << endl;
  cerr << "     2. position         "    << endl;
  cerr << "     3. pFst probability "    << endl  << endl;

  cerr << "INFO: usage:  pFst --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf --deltaaf 0.1 --type PL" << endl;
  cerr << endl;
  cerr << "INFO: required: t,target     -- argument: a zero based comma separated list of target individuals corresponding to VCF columns       "  << endl;
  cerr << "INFO: required: b,background -- argument: a zero based comma separated list of background individuals corresponding to VCF columns   "  << endl;
  cerr << "INFO: required: f,file       -- argument: a properly formatted VCF.                                                                  "  << endl;
  cerr << "INFO: required: y,type       -- argument: genotype likelihood format ; genotypes: GP, GL or PL; pooled: PO                           "  << endl;
  cerr << "INFO: optional: d,deltaaf    -- argument: skip sites where the difference in allele frequencies is less than deltaaf, default is zero"  << endl;
  cerr << "INFO: optional: r,region     -- argument: a tabix compliant genomic range : seqid or seqid:start-end                                 "  << endl;
  cerr << "INFO: optional: c,counts     -- switch  : use genotype counts rather than genotype likelihoods to estimate parameters, default false "  << endl;

  cerr << endl;

  printVersion() ;
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

double logLbinomial(double x, double n, double p){

  double ans = lgamma(n+1)-lgamma(x+1)-lgamma(n-x+1) + x * log(p) + (n-x) * log(1-p);
  return ans;
    
}

int main(int argc, char** argv) {

  // pooled or genotyped

  int pool = 0;

  // the filename

  string filename = "NA";

  // set region to scaffold

  string region = "NA"; 

  // using vcflib; thanks to Erik Garrison 

  VariantCallFile variantFile;

  // zero based index for the target and background indivudals 
  
  map<int, int> it, ib;
  
  // deltaaf is the difference of allele frequency we bother to look at 

  string deltaaf ;
  double daf  = 0;

  // 

  int counts = 0;

  // type pooled GL PL

  string type = "NA";

    const struct option longopts[] = 
      {
	{"version"   , 0, 0, 'v'},
	{"help"      , 0, 0, 'h'},
	{"counts"    , 0, 0, 'c'},
        {"file"      , 1, 0, 'f'},
	{"target"    , 1, 0, 't'},
	{"background", 1, 0, 'b'},
	{"deltaaf"   , 1, 0, 'd'},
	{"region"    , 1, 0, 'r'},
	{"type"      , 1, 0, 'y'},
	{0,0,0,0}
      };

    int index;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "r:d:t:b:f:y:chv", longopts, &index);
	
	switch (iarg)
	  {
	  case 'h':
	    printHelp();
	    return 0;
	  case 'v':
	    printVersion();
	    return 0;
	  case 'y':	    
	    type = optarg;
	    cerr << "INFO: genotype likelihoods set to: " << type << endl;
	    if(type == "GT"){
	      cerr << "INFO: using counts flag as GT was specified" << endl;
	    }
	    break;
	  case 'c':
	    cerr << "INFO: using genotype counts rather than genotype likelihoods" << endl;
	    counts = 1;
	    break;
	  case 't':
	    loadIndices(ib, optarg);
	    cerr << "INFO: There are " << ib.size() << " individuals in the target" << endl;
	    cerr << "INFO: target ids: " << optarg << endl;
	    break;
	  case 'b':
	    loadIndices(it, optarg);
	    cerr << "INFO: There are " << it.size() << " individuals in the background" << endl;
	    cerr << "INFO: background ids: " << optarg << endl;
	    break;
	  case 'f':
	    cerr << "INFO: File: " << optarg  <<  endl;
	    filename = optarg;
	    break;
	  case 'd':
	    cerr << "INFO: only scoring sites where the allele frequency difference is greater than: " << optarg << endl;
	    deltaaf = optarg;
	    daf = atof(deltaaf.c_str());	    
	    break;
	  case 'r':
            cerr << "INFO: set seqid region to : " << optarg << endl;
	    region = optarg; 
	    break;
	  default:
	    break;
	  }
      }

    
    if(filename == "NA"){
      cerr << "FATAL: did not specify the file\n";
      printHelp();
      exit(1);
    }
    
    variantFile.open(filename);
    
    if(region != "NA"){
      if(! variantFile.setRegion(region)){
	cerr <<"FATAL: unable to set region" << endl;
	return 1;
      }
    }

    if (!variantFile.is_open()) {
      exit(1);
    }
    map<string, int> okayGenotypeLikelihoods;
    okayGenotypeLikelihoods["PL"] = 1;
    okayGenotypeLikelihoods["GT"] = 1;
    okayGenotypeLikelihoods["PO"] = 1;
    okayGenotypeLikelihoods["GL"] = 1;
    okayGenotypeLikelihoods["GP"] = 1;

    if(type == "GT"){
      counts = 1;
    }

    if(type == "NA"){
      cerr << "FATAL: failed to specify genotype likelihood format : PL,PO,GL,GP" << endl;
      printHelp();
      return 1;
    }
    if(okayGenotypeLikelihoods.find(type) == okayGenotypeLikelihoods.end()){
      cerr << "FATAL: genotype likelihood is incorrectly formatted, only use: PL,PO,GL,GP" << endl;
      printHelp();
      return 1;
    }    

    Variant var(variantFile);

    vector<string> samples = variantFile.sampleNames;
    int nsamples = samples.size();

    while (variantFile.getNextVariant(var)) {

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
		total.push_back(sample);		
	      }
	      if(ib.find(index) != ib.end()){
		background.push_back(sample);
		total.push_back(sample);		
	      }
	    }
            
	index += 1;
	}
	
	zvar * populationTarget        ;
	zvar * populationBackground    ;
	zvar * populationTotal         ;

	if(type == "PO"){
	  populationTarget     = new pooled();
          populationBackground = new pooled();
          populationTotal      = new pooled();
	}
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
	if(type == "GT"){
          populationTarget     = new gt();
          populationBackground = new gt();
          populationTotal      = new gt();
        }	

	populationTotal->loadPop(total          , var.sequenceName, var.position);	
	populationTarget->loadPop(target        , var.sequenceName, var.position);
	populationBackground->loadPop(background, var.sequenceName, var.position);

	if(populationTarget->npop < 2 || populationBackground->npop < 2){
          delete populationTarget;
	  delete populationBackground;
          delete populationTotal;

	  continue;
	}

	populationTotal->estimatePosterior();	
	populationTarget->estimatePosterior();
	populationBackground->estimatePosterior();

	if(populationTarget->alpha == -1 || populationBackground->alpha == -1){
	  delete populationTarget;
	  delete populationBackground;
	  delete populationTotal;


          continue;
        }

	if(counts == 1){

	  populationTotal->alpha  = 0.001 + populationTotal->nref;
	  populationTotal->beta   = 0.001 + populationTotal->nalt;
	  
	  populationTarget->alpha = 0.001 + populationTarget->nref;
	  populationTarget->beta  = 0.001 + populationTarget->nalt;

	  populationBackground->alpha = 0.001 + populationBackground->nref;
	  populationBackground->beta  = 0.001 + populationBackground->nalt;


	}

	double populationTotalEstAF       = bound(populationTotal->beta      / (populationTotal->alpha      + populationTotal->beta)     );
	double populationTargetEstAF      = bound(populationTarget->beta     / (populationTarget->alpha     + populationTarget->beta)    );
	double populationBackgroundEstAF  = bound(populationBackground->beta / (populationBackground->alpha + populationBackground->beta));

	// cout << populationTotalEstAF << "\t" << populationTotal->af << endl;

	// x, n, p
	double null = logLbinomial(populationTarget->beta, (populationTarget->alpha + populationTarget->beta),  populationTotalEstAF) +
	  logLbinomial(populationBackground->beta, (populationBackground->alpha + populationBackground->beta),  populationTotalEstAF) ;
	double alt  = logLbinomial(populationTarget->beta, (populationTarget->alpha + populationTarget->beta),  populationTargetEstAF) +
	  logLbinomial(populationBackground->beta, (populationBackground->alpha + populationBackground->beta),  populationBackgroundEstAF) ;
		
	double l = 2 * (alt - null);
	
	if(l <= 0){
	  delete populationTarget;
	  delete populationBackground;
	  delete populationTotal;


	  continue;
	}

	int     which = 1;
	double  p ;
	double  q ;
	double  x  = l;
	double  df = 1;
	int     status;
	double  bound ;
	
	cdfchi(&which, &p, &q, &x, &df, &status, &bound );
	
	cout << var.sequenceName << "\t"  << var.position << "\t" << 1-p << endl ;
	
	delete populationTarget;
	delete populationBackground;
	delete populationTotal;
	
	populationTarget     = NULL;
	populationBackground = NULL;
	populationTotal      = NULL;

    }
    return 0;		    
}
