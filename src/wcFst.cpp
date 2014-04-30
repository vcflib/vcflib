#include "Variant.h"
#include "split.h"
#include "cdflib.hpp"
#include "pdflib.hpp"

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

struct pop{

  double nalt ;
  double nref ;
  double af   ; 
  double nhomr;
  double nhoma;
  double nhet ;
  double ngeno;
  double hfrq ; 
  double fis  ;

  vector<int>    geno_index ;
  vector< vector< double > > unphred_p;  
  
};

double unphred(string phred){
  
  double unphred = atof(phred.c_str());
  unphred = unphred / -10;
  return unphred;
  
}

void initPop(pop & population){
  population.nalt  = 0;
  population.nref  = 0;
  population.af    = 0;
  population.nhomr = 0;
  population.nhoma = 0;
  population.nhet  = 0;
  population.ngeno = 0;
  population.fis   = 0;
  population.hfrq  = 0; 
}

void loadPop( vector< map< string, vector<string> > >& group, pop & population){

  vector< map< string, vector<string> > >::iterator targ_it = group.begin();
  
  for(; targ_it != group.end(); targ_it++){
    
    
    string genotype = (*targ_it)["GT"].front();
    
    while(1){
      if(genotype == "0/0"){
	population.ngeno += 1;
	population.nhomr += 1;
	population.nref  += 2;
	population.geno_index.push_back(0);
	break;
      }
      if(genotype == "0/1"){
	population.ngeno += 1;
	population.nhet  += 1;
	population.nref  += 1;
	population.nalt  += 1;
	population.geno_index.push_back(1);
	break;
      }
      if(genotype == "1/1"){
	population.ngeno += 1;
	population.nhoma += 1;
	population.nalt  += 2;
	population.geno_index.push_back(2);
	break;
      }
      if(genotype == "0|0"){
	population.ngeno += 1;
	population.nhomr += 1;
	population.nref  += 2;
	population.geno_index.push_back(0);
	break;
      }
      if(genotype == "0|1"){
	population.ngeno += 1;
	population.nhet  += 1;
	population.nref  += 1;
	population.nalt  += 1;
	population.geno_index.push_back(1);
	break;
      }
      if(genotype == "1|0"){
        population.ngeno += 1;
        population.nhet  += 1;
        population.nref  += 1;
        population.nalt  += 1;
	population.geno_index.push_back(1);
        break;
      }
      if(genotype == "1|1"){
	population.ngeno += 1;
	population.nhoma += 1;
	population.nalt  += 2;
	population.geno_index.push_back(2);
	break;
      }
      if(genotype == "1/0"){
        population.ngeno += 1;
        population.nhet  += 1;
        population.nref  += 1;
        population.nalt  += 1;
        population.geno_index.push_back(1);
        break;
      }
      cerr << "FATAL: unknown genotype:" <<  genotype << endl;
      exit(1);
    }
  }
  
  if(population.nalt == 0 && population.nref == 0){
    population.af = -1;
  }
  else{
    
    population.af  = (population.nalt / (population.nref + population.nalt));
    population.hfrq = population.nhet / population.ngeno;
    
    if(population.nhet > 0){
      population.fis = ( 1 - ((population.nhet/population.ngeno) / (2*population.af*(1 - population.af))));
    }
    else{
      population.fis = 1;
    }
    if(population.fis < 0){
      population.fis = 0.00001;
    }
  }  
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

void getPosterior(pop & population, double *alpha, double *beta){
  
  int ng = population.geno_index.size();

  for(int i = 0 ; i < ng; i++){
    
    double aa = population.unphred_p[i][0] ; 
    double ab = population.unphred_p[i][1] ; 
    double bb = population.unphred_p[i][2] ; 

    double norm = log(exp(aa) + exp(ab) + exp(bb));

    int gi = population.geno_index[i];
    
    (*alpha) += exp(ab - norm);
    (*beta ) += exp(ab - norm);

    (*alpha) += 2 * exp(aa - norm);
    (*beta ) += 2 * exp(bb - norm);
    
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

  string deltaaf ;
  double daf  = 0;

  // 


    const struct option longopts[] = 
      {
	{"version"   , 0, 0, 'v'},
	{"help"      , 0, 0, 'h'},
        {"file"      , 1, 0, 'f'},
	{"target"    , 1, 0, 't'},
	{"background", 1, 0, 'b'},
	{"deltaaf"   , 1, 0, 'd'},
	{"region"    , 1, 0, 'r'},
	{0,0,0,0}
      };

    int index;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "r:d:t:b:f:chv", longopts, &index);
	
	switch (iarg)
	  {
	  case 'h':
	    cerr << endl << endl;
	    cerr << "INFO: help" << endl;
	    cerr << "INFO: description:" << endl;
	    cerr << "      wcFst is Weir & Cockerham's Fst for two populations.  Negative values are VALID,  " << endl;
	    cerr << "      they are sites which can be treated as zero Fst. For more information see Evolution, Vol. 38 N. 6 Nov 1984. " << endl;
	    cerr << "      Specifically wcFst uses equations 1,2,3,4.                                                              " << endl << endl;

	    cerr << "Output : 3 columns :     "    << endl;
	    cerr << "     1. seqid            "    << endl;
	    cerr << "     2. position         "    << endl;
	    cerr << "     3. wcFst            "    << endl  << endl;

	    cerr << "INFO: usage:  wcFst --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf --deltaaf 0.1" << endl;
	    cerr << endl;
	    cerr << "INFO: required: t,target     -- a zero based comma seperated list of target individuals corrisponding to VCF columns" << endl;
	    cerr << "INFO: required: b,background -- a zero based comma seperated list of background individuals corrisponding to VCF columns" << endl;
	    cerr << "INFO: required: f,file a     -- proper formatted VCF.  the FORMAT field MUST contain \"PL\"" << endl; 
	    cerr << "INFO: optional: d,deltaaf    -- skip sites where the difference in allele frequencies is less than deltaaf, default is zero"      << endl;
	    cerr << endl; 
	    cerr << "INFO: version 1.0.0 ; date: April 2014 ; author: Zev Kronenberg; email : zev.kronenberg@utah.edu " << endl;
	    cerr << endl << endl;
	    return 0;
	  case 'v':
	    cerr << endl << endl;
	    cerr << "INFO: version 1.0.0 ; date: April 2014 ; author: Zev Kronenberg; email : zev.kronenberg@utah.edu "  << endl;
	    return 0;
	  case 't':
	    loadIndices(it, optarg);
	    cerr << "INFO: There are " << it.size() << " individuals in the target" << endl;
	    cerr << "INFO: target ids: " << optarg << endl;
	    break;
	  case 'b':
	    loadIndices(ib, optarg);
	    cerr << "INFO: There are " << ib.size() << " individuals in the background" << endl;
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
    
    variantFile.open(filename);
    
    if(region != "NA"){
      variantFile.setRegion(region);
    }

    if (!variantFile.is_open()) {
        return 1;
    }
    
    Variant var(variantFile);

    while (variantFile.getNextVariant(var)) {
        map<string, map<string, vector<string> > >::iterator s     = var.samples.begin(); 
        map<string, map<string, vector<string> > >::iterator sEnd  = var.samples.end();
        
	// biallelic sites naturally 

	if(var.alt.size() > 1){
	  continue;
	}
	
	vector < map< string, vector<string> > > target, background, total;
	        
	int index = 0;

        for (; s != sEnd; ++s) {

            map<string, vector<string> >& sample = s->second;

	    if(sample["GT"].front() != "./."){
	      if(it.find(index) != it.end() ){
		target.push_back(sample);
	      }
	      if(ib.find(index) != ib.end()){
		background.push_back(sample);
	      }
	    }            
	    index += 1;
	}
	
	if(target.size() < 5 || background.size() < 5 ){
	  continue;
	}
	
	pop popt, popb;
	
	initPop(popt);
	initPop(popb);

	loadPop(target,     popt);
	loadPop(background, popb);


	if(popt.af == -1 || popb.af == -1){
	  continue;
	}


	double afdiff = abs(popt.af - popb.af);

	if(afdiff < daf){
	  continue;
	}

	
	// pg 1360 B.S Weir and C.C. Cockerham 1984
	double nbar = ( popt.ngeno / 2 ) + (popb.ngeno / 2);
	double rn   = 2*nbar;
	
	// special case of only two populations
	double nc   =  rn ;
	nc -= (pow(popt.ngeno,2)/rn);
	nc -= (pow(popb.ngeno,2)/rn);
	// average sample frequency
	double pbar = (popt.af + popb.af) / 2;

	//	pbar += ((popt.af * popt.ngeno)/ rn);
	//pbar += ((popb.af * popb.ngeno)/ rn);
	// sample variance of allele A frequences over the population 
	
	double s2 = 0;
	s2 += ((popt.ngeno * pow(popt.af - pbar, 2))/nbar);
	s2 += ((popb.ngeno * pow(popb.af - pbar, 2))/nbar);
	
	// average heterozygosity 
	
	double hbar = (popt.hfrq + popb.hfrq) / 2;
	//	hbar += ((popt.ngeno*popt.hfrq)/rn);
	//      hbar += ((popb.ngeno*popb.hfrq)/rn);
	
	//global af var
	double pvar = pbar * (1 - pbar);

	// a, b, c

	double avar1 = nbar / nc;
	double avar2 = 1 / (nbar -1) ;
	double avar3 = pvar - (0.5*s2) - (0.25*hbar);
	double avar  = avar1 * (s2 - (avar2 * avar3));
	
	double bvar1 = nbar / (nbar - 1);
	double bvar2 = pvar - (0.5*s2) - (((2*nbar -1)/(4*nbar))*hbar);
	double bvar  = bvar1 * bvar2;

	double cvar = 0.5*hbar;
	
	double fst = avar / (avar+bvar+cvar);
	
	cout << var.sequenceName << "\t"  << var.position << "\t" << popt.af << "\t" << popb.af << "\t" << fst << endl ;

    }
    return 0;		    
}
