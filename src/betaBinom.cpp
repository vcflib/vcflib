#include "Variant.h"
#include "split.h"
#include "pdflib.hpp"
#include "cdflib.hpp"

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
  
}

void loadPop( vector< map< string, vector<string> > >& group, pop & population){

  vector< map< string, vector<string> > >::iterator targ_it = group.begin();
  
  for(; targ_it != group.end(); targ_it++){
    
    population.ngeno += 1;
    
    string genotype = (*targ_it)["GT"].front();
    
    vector<double> phreds;

	phreds.push_back( unphred((*targ_it)["PL"][0]));
	phreds.push_back( unphred((*targ_it)["PL"][1]));
	phreds.push_back( unphred((*targ_it)["PL"][2]));
    
	population.unphred_p.push_back(phreds);
    
	while(1){
	  if(genotype == "0/0"){
	    population.nhomr += 1;
	    population.nref  += 2;
	    population.geno_index.push_back(0);
	    break;
	  }
	  if(genotype == "0/1"){
	    population.nhet  += 1;
	    population.nref  += 1;
	    population.nalt  += 1;
	    population.geno_index.push_back(1);
	    break;
	  }
	  if(genotype == "1/1"){
	    population.nhoma += 1;
	    population.nalt  += 2;
	    population.geno_index.push_back(2);
	    break;
	  }
	  if(genotype == "0|0"){
	    population.nhomr += 1;
	    population.nref  += 2;
	    population.geno_index.push_back(0);
	    break;
	  }
	  if(genotype == "0|1"){
	    population.nhet  += 1;
	    population.nref  += 1;
	    population.nalt  += 1;
	    population.geno_index.push_back(1);
	    break;
	  }
	  if(genotype == "1|1"){
	    population.nhoma += 1;
	    population.nalt  += 2;
	    population.geno_index.push_back(2);
	    break;
	  }
	  cerr << "FATAL: unknown genotype" << endl;
	  exit(1);
	}
  }

  if(population.nalt == 0 && population.nref == 0){
    population.af = -1;
  }
  else{
   
    population.af  = (population.nalt / (population.nref + population.nalt));

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

  double counts = 1;

    const struct option longopts[] = 
      {
	{"version"   , 0, 0, 'v'},
	{"help"      , 0, 0, 'h'},
        {"file"      , 1, 0, 'f'},
	{"target"    , 1, 0, 't'},
	{"background", 1, 0, 'b'},
	{"deltaaf"   , 1, 0, 'd'},
	{"region"    , 1, 0, 'r'},
	{"genotype"  , 1, 0, 'g'},
	{0,0,0,0}
      };

    int index;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "r:d:t:b:f:hv", longopts, &index);
	
	switch (iarg)
	  {
	  case 'h':
	    cerr << endl << endl;
	    cerr << "INFO: help" << endl;
	    cerr << "INFO: description: bFst is a bayesian approach to Fst, where Fst is a free parameter.  Importantly bFst taken genotype  " << endl;
	    cerr << "                   into account, meaning sites of poor quality will result in wider credible intervals.  For more more  " << endl;
	    cerr << "                   information see \" A bayesian approach to inferring population structure from dominant markers. bFst " << endl;
	    cerr << "                   is not using the likelihood in equation 3.                                                           " << endl;

	    cerr << "INFO: usage:  bFst --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf --deltaaf 0.1" << endl;
	    cerr << endl;
	    cerr << "INFO: required: t,target     -- a zero bases comma seperated list of target individuals corrisponding to VCF columns" << endl;
	    cerr << "INFO: required: b,background -- a zero bases comma seperated list of target individuals corrisponding to VCF columns" << endl;
	    cerr << "INFO: required: f,file a     -- proper formatted VCF.  the FORMAT field MUST contain \"PL\"" << endl; 
	    cerr << "INFO: required: d,deltaaf    -- skip sites were the difference in allele frequency is less than deltaaf" << endl;
	    cerr << endl; 
	    cerr << "INFO: version 1.0.0 ; date: April 2014 ; author: Zev Kronenberg; email : zev.kronenberg@utah.edu " << endl;
	    cerr << endl << endl;
	    return 0;

	  case 'v':
	    cerr << endl << endl;
	    cerr << "INFO: version 1.0.0 ; date: April 2014 ; author: Zev Kronenberg; email : zev.kronenberg@utah.edu "  << endl;
	    cerr << endl << endl;
	    return 0;

	  case 't':
	    loadIndices(ib, optarg);
	    cerr << endl << endl;
	    cerr << "INFO: There are " << ib.size() << " individuals in the target" << endl;
	    cerr << endl << endl;
	    break;

	  case 'b':
	    loadIndices(it, optarg);
	    cerr << endl << endl;
	    cerr << "INFO: There are " << it.size() << " individuals in the background" << endl;
	    cerr << endl << endl;
	    break;

	  case 'f':
	    cerr << endl << endl;
	    cerr << "INFO: File: " << optarg  <<  endl;
	    cerr << endl << endl;
	    filename = optarg;
	    break;
	  case 'd':
	    cerr << endl << endl;
	    cerr << "INFO: difference in allele frequency : " << optarg << endl;
	    cerr << endl << endl;
	    deltaaf = optarg;
	    daf = atof(deltaaf.c_str());	    
	    break;
	  case 'r':
	    cerr << endl << endl;
            cerr << "INFO: set seqid region to : " << optarg << endl;
	    region = optarg; 
	    break;
	  default:
	    printf ("?? getopt returned character code 0%o ??\n", iarg);
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
		total.push_back(sample);
		
	      }
	      if(ib.find(index) != ib.end()){
		background.push_back(sample);
		total.push_back(sample);
	      }
	    }
            
	index += 1;
	}
	
	if(target.size() < 5 || background.size() < 5 ){
	  continue;
	}
	
	pop popt, popb, popTotal;
	
	initPop(popt);
	initPop(popb);
	initPop(popTotal);

	loadPop(target,     popt);
	loadPop(background, popb);
	loadPop(total,  popTotal);

	if(popt.af == -1 || popb.af == -1){
	  continue;
	}
	if(popt.af == 1  && popb.af == 1){
	  continue;
	}
	if(popt.af == 0 && popb.af  == 0){
	  continue;
	}

	double afdiff = abs(popt.af - popb.af);

	if(afdiff < daf){
	  continue;
	}

	double alphaT = 0.01;
	double alphaB = 0.01;
	double betaT  = 0.01;
	double betaB  = 0.01;
      
	getPosterior(popt, &alphaT, &betaT);
	getPosterior(popb, &alphaB, &betaB);
	
	double targm = alphaT / ( alphaT + betaT );
	double backm = alphaB / ( alphaB + betaB );
	
    
	double p, q, x, y, bound;
	p = -1;
	q = -1;
	
	int status, which;
	which = 1;

	x = targm;
	y = 1 - x; 

	cdfbet(&which, &p, &q, &x, &y, &alphaB, &betaB, &status, &bound);



	int n = int(popt.nalt + popt.nref);

	
	cout << var.sequenceName << "\t"  << var.position << "\t" << p << endl ;

    }
    return 0;		    
}
