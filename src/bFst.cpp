#include "Variant.h"
#include "split.h"
#include "pdflib.hpp"

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

struct pop{

  double nalt ;
  double nref ;
  double af   ; 
  double nhomr;
  double nhoma;
  double nhet ;
  double ngeno;
  double fis  ;
  vector<int>    questionable;
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

  int index = 0;

  for(; targ_it != group.end(); targ_it++){
    
    string genotype = (*targ_it)["GT"].front();
    
    vector<double> phreds;

    phreds.push_back( unphred((*targ_it)["PL"][0]));
    phreds.push_back( unphred((*targ_it)["PL"][1]));
    phreds.push_back( unphred((*targ_it)["PL"][2]));
    
    double scaled ;
    double norm   = log(exp(phreds[0]) + exp(phreds[1]) + exp(phreds[2]));
    
    population.unphred_p.push_back(phreds);
    
    while(1){
      if(genotype == "0/0"){
	population.ngeno += 1;
	population.nhomr += 1;
	population.nref  += 2;
	population.geno_index.push_back(0);	    
	scaled = exp(phreds[0] - norm); 
	break;
      }
      if(genotype == "0/1"){
	population.ngeno += 1;
	population.nhet  += 1;
	population.nref  += 1;
	population.nalt  += 1;
	population.geno_index.push_back(1);
	scaled = exp(phreds[1] - norm); 
	break;
      }
      if(genotype == "1/1"){
	population.ngeno += 1;
	population.nhoma += 1;
	population.nalt  += 2;
	population.geno_index.push_back(2);
	scaled = exp(phreds[2] - norm); 
	break;
      }
      if(genotype == "0|0"){
	population.ngeno += 1;
	population.nhomr += 1;
	population.nref  += 2;
	population.geno_index.push_back(0);
	scaled = exp(phreds[0] - norm); 
	break;
      }
      if(genotype == "0|1"){
	population.ngeno += 1;
	population.nhet  += 1;
	population.nref  += 1;
	population.nalt  += 1;
	population.geno_index.push_back(1);
	scaled = exp(phreds[1] - norm); 
	break;
      }
      if(genotype == "1|1"){
	population.ngeno += 1;
	population.nhoma += 1;
	population.nalt  += 2;
	population.geno_index.push_back(2);
	scaled = exp(phreds[2] - norm); 
	break;
      }
      cerr << "FATAL: unknown genotype" << endl;
      exit(1);
    }
    
    if(scaled < 0.75){
      population.questionable.push_back(index);
    }	
    index += 1;
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


void phardy(vector<double>& results, double af, double fis){
  
  double p0 = pow((1 - af),2) + ((1 - af)*af*fis);
  double p1 = 2*(1 - af)*af*(1 - fis);
  double p2 = pow(af,2) + ((1 - af)*af*fis);
  
  results.push_back(p0);
  results.push_back(p1);
  results.push_back(p2);
  
}

double likelihood(pop & population, double af, double fis){

  af  = bound(af);
  fis = bound(fis);

  double het;
  double all;
  double ref;
  double alt;
  double loglikelihood = 0;

  vector<double> genotypeProbs;

  phardy(genotypeProbs, af, fis);

  vector<int>::iterator it = population.geno_index.begin();

  int geno_indx = 0;

  for(; it != population.geno_index.end(); it++){

    double aa = population.unphred_p[geno_indx][0] + log(genotypeProbs[0]);
    double ab = population.unphred_p[geno_indx][1] + log(genotypeProbs[1]);
    double bb = population.unphred_p[geno_indx][2] + log(genotypeProbs[2]);
    double norm = exp(aa) + exp(ab) + exp(bb);

    double prop = population.unphred_p[geno_indx][*it] +  log(genotypeProbs[*it]);
    
    loglikelihood += (prop - norm);

    geno_indx++;
  }
  
  return loglikelihood;
  
}

double FullProb(pop & target, pop & back, vector<double>& p)
{

// parameters targetAf backgroundAf targetFis backgroundFis totalAf fst

  double alpha = ( (1-p[5])/p[5] ) * p[4];
  double beta  = ( (1-p[5])/p[5] ) * (1 - p[4]);


  double afprior  = log( r8_normal_pdf (p[6], 0.1, p[4]));
  double afpriorT = log( r8_normal_pdf (target.af, 0.05, p[0]));
  double afpriorB = log( r8_normal_pdf (back.af,   0.05, p[1]));
  
  if(std::isinf(afprior) || std::isnan(afprior)){
    return -100000;
  }

  double ptaf = log( r8_beta_pdf(alpha, beta, p[0]) );
  double pbaf = log( r8_beta_pdf(alpha, beta, p[1]) );

  if( std::isinf(ptaf) || std::isnan(ptaf) || std::isinf(pbaf) || std::isnan(pbaf) ){
    return -100000;
  }


  double llt  = likelihood(target, p[0], p[2]);
  double llb  = likelihood(back,   p[1], p[3]);
  double full = llt + llb + ptaf + pbaf + afprior + afpriorT + afpriorB;

  
  return full;
  
}



void updateParameters(pop & target, pop & background, vector<double>& parameters, int pindx){

  // parameters targetAf backgroundAf targetFis backgroundFis totalAf fst

  double origpar = parameters[pindx];

  double accept  = ((double)rand() / (double)(RAND_MAX));
  double up      = ((double)rand() / (double)(RAND_MAX))/10 - 0.05;
  double updatep = parameters[pindx] + up;
   
  //  cerr << accept << "\t" << up << endl;

  if(updatep >= 1 || updatep <= 0){
    return;
  }
 
  double llB = FullProb(target, background, parameters);
  parameters[pindx] = updatep;
  double llT = FullProb(target, background, parameters);
  
  if((llT - llB) > accept){
    return; 
  }
  else{
    parameters[pindx] = origpar;
  }
} 

void updateGenotypes(pop & target, pop & background, vector<double>& parameters, int gindex, int tbindex){
  
  // tbindex indicates if the subroutine will update the target or background genotype;

  double accept  = ((double)rand() / (double)(RAND_MAX));
  int  newGindex = rand() % 3;
  
  //cerr << newGindex << endl;
  //cerr << "gindex "   << gindex << endl;
  //cerr << "gsize t:" << target.geno_index.size() << endl;  
  //cerr << "gsize b:" << background.geno_index.size() << endl;  
  


  int oldtindex = target.geno_index[gindex] ;
  int oldbindex = background.geno_index[gindex] ;  
  
  double llB = FullProb(target, background, parameters);
  
  if(tbindex == 0){
    //udate target                                                                                                  
    target.geno_index[gindex] = newGindex;
      }
  else{
    // update background
    background.geno_index[gindex] = newGindex;
  }
  double llT = FullProb(target, background, parameters);
  
  if((llT - llB) > accept){
    return;
  }
  else{
    if(tbindex == 0){
      target.geno_index[gindex] = oldtindex;
    }
    else{
      target.geno_index[gindex] = oldbindex;
    }
  } 
}


int cmp(const void *x, const void *y)
{
  double xx = *(double*)x, yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
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

  // using vcflib; thanks to Erik Garrison 
  
  VariantCallFile variantFile ;

  // zero based index for the target and background indivudals 
  
  map<int, int> it, ib;
  
  // deltaaf is the difference of allele frequency we bother to look at 

  string deltaaf ;
  double daf  = -1;

    const struct option longopts[] = 
      {
	{"version"   , 0, 0, 'v'},
	{"help"      , 0, 0, 'h'},
        {"file"      , 1, 0, 'f'},
	{"target"    , 1, 0, 't'},
	{"background", 1, 0, 'b'},
	{"deltaaf"   , 1, 0, 'd'},
	{0,0,0,0}
      };

    int index;
    int iarg = 0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "d:t:b:f:hv", longopts, &index);
	
	switch (iarg)
	  {
	  case 0:
	    break;
	  case 'h':
	    cerr << endl;

	    cerr << "INFO: help: " << endl << endl;

	    cerr << "     bFst is a Bayesian approach to Fst.  Importantly bFst account for genotype uncertainty in the model using genotype likelihoods."       << endl;
	    cerr << "     For a more detailed description see: Holsinger et al. Molecular Ecology Vol 11, issue 7 2002.  The likelihood function has been "	     << endl;
	    cerr << "     modified to use genotype likelihoods provided by variant callers. There are five free parameters estimated in the model: each "	     << endl;
	    cerr << "     subpopulation's allele frequency and Fis (fixation index, within each subpopulation), a free parameter for the total population\'s "  << endl;
	    cerr << "     allele frequency, and Fst. "                                                                                      << endl             << endl;
	
	      cerr << "Output : 11 columns :                          " << endl; 
	      cerr << "     1.  Seqid                                     " << endl;
	      cerr << "     2.  Position				     " << endl;
	      cerr << "     3.  Observed allele frequency in target.	     " << endl;
	      cerr << "     4.  Estimated allele frequency in target.     " << endl;
	      cerr << "     5.  Observed allele frequency in background.  " << endl;
	      cerr << "     6.  Estimated allele frequency in background. " << endl;
	      cerr << "     7.  Observed allele frequency combined. 	     " << endl;
	      cerr << "     8.  Estimated allele frequency in combined.   " << endl;
	      cerr << "     9.  ML estimate of Fst (mean)		     " << endl;
	      cerr << "     10. Lower bound of the 95% credible interval  " << endl;
	      cerr << "     11. Upper bound of the 95% credible interval  " << endl << endl;
											 

	    cerr << "INFO: usage:  bFst --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf --deltaaf 0.1" << endl;
	    cerr << endl;
	    cerr << "INFO: required: t,target     -- a zero bases comma separated list of target individuals corrisponding to VCF columns" << endl;
	    cerr << "INFO: required: b,background -- a zero bases comma separated list of background individuals corrisponding to VCF columns" << endl;
	    cerr << "INFO: required: f,file a     -- a proper formatted VCF file.  the FORMAT field MUST contain \"PL\"" << endl; 
	    cerr << "INFO: required: d,deltaaf    -- skip sites were the difference in allele frequency is less than deltaaf" << endl;
	    cerr << endl; 
	    printVersion();
	    cerr << endl << endl;
	    return 0;

	  case 'v':
	    printVersion();
	    return 0;

	  case 't':
	    loadIndices(ib, optarg);
	    cerr << "INFO: There are " << ib.size() << " individuals in the target" << endl;
	    break;

	  case 'b':
	    loadIndices(it, optarg);
	    cerr << "INFO: There are " << it.size() << " individuals in the background" << endl;
	    break;

	  case 'f':
	    cerr << "INFO: File: " << optarg  <<  endl;
	    filename = optarg;
	    break;

	  case 'd':
	    cerr << "INFO: difference in allele frequency : " << optarg << endl;
	    deltaaf = optarg;
	    daf = atof(deltaaf.c_str());	    
	    break;
	  default: 
	    break; 
	    cerr << endl;
	    cerr << "FATAL: unknown command line option " << optarg << endl << endl ;
	    cerr << "INFO:  please use bFst --help      " << endl; 
	    cerr << endl;
	    return(1);
	  }

      }

    if(daf == -1){
    cerr << endl;
      cerr << "FATAL: did not specify deltaaf" << endl;
      cerr << "INFO:  please use bFst --help      " << endl; 
      cerr << endl;
      return(1);
    }

    if(filename == "NA"){
      cerr << endl;
      cerr << "FATAL: did not specify VCF file" << endl;
      cerr << "INFO:  please use bFst --help      " << endl; 
      cerr << endl;
      return(1);
    }

    variantFile.open(filename);
    

    if (!variantFile.is_open()) {
      cerr << endl;
      cerr << "FATAL: could not open VCF file" << endl;
      cerr << "INFO:  please use bFst --help" << endl; 
      cerr << endl;
      return(1);
    }
    if(it.size() < 2){
      cerr << endl;
      cerr << "FATAL: target not specified or less than two indviduals" << endl; 
      cerr << "INFO:  please use bFst --help                          " << endl; 
      cerr << endl;
    }
    if(ib.size() < 2){
      cerr << endl;
      cerr << "FATAL: target not specified or less than two indviduals"<< endl;
      cerr << "INFO:  please use bFst --help                          " << endl;
      cerr << endl;
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
	      total.push_back(sample);
	      
	    }
	    if(ib.find(index) != ib.end()){
		background.push_back(sample);
		total.push_back(sample);
	    }
	  }
    
	  index += 1;
	}
	
	if(target.size() < 2 || background.size() < 2 ){
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
	
	
	cerr << "INFO: target has "     << popt.questionable.size() << " questionable genotypes " << endl;
	cerr << "INFO: background has " << popb.questionable.size() << " questionable genotypes " << endl;

	// Parameters- targetAf backgroundAf targetFis backgroundFis totalAf fst
	vector<double> parameters;
	parameters.push_back(popt.af);
	parameters.push_back(popb.af);
	parameters.push_back(popt.fis);
	parameters.push_back(popb.fis);
	parameters.push_back(popTotal.af);
	parameters.push_back(0.1);
	parameters.push_back(popTotal.af);

	double sums [6] = {0};
	double fsts [10000]  ;

	for(int i = 0; i < 15000; i++){
	  
	  // update each of j parameters
	  
	  for(int j = 0; j < 6; j++ ){
	    
	    updateParameters(popt, popb, parameters, j);
	    if(i > 4999){
	      sums[j]     += parameters[j]; 
	    }
	  }
	  if(i > 4999){
	    fsts[i - 5000] =  parameters[5]; 
	  }
	  for(vector<int>::iterator itt = popt.questionable.begin(); itt != popt.questionable.end(); itt++){
	    updateGenotypes(popt, popb, parameters, (*itt), 0);

	  }
	  for(vector<int>::iterator itb = popb.questionable.begin(); itb != popb.questionable.end(); itb++){
	    updateGenotypes(popt, popb, parameters, (*itb) , 1);
	  }
	}
		
	qsort (fsts, sizeof(fsts)/sizeof(fsts[0]), sizeof(fsts[0]), cmp );
	
	double lcredint = fsts[500];
	double hcredint = fsts[9500]; 
	
    	cout << var.sequenceName << "\t"  << var.position     
	     << "\t"  << popt.af
             << "\t"  << sums[0]/10000
	     << "\t"  << popb.af 
	     << "\t"  << sums[1]/10000
	     << "\t"  << popTotal.af 
	     << "\t"  << sums[4]/10000
	     << "\t"  << sums[5]/10000
	     << "\t"  << lcredint
	     << "\t"  << hcredint
	     << endl;
    }
    return 0;		    
}
