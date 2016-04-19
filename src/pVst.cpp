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
#include <algorithm> 
#include <ctime>        
#include <cstdlib>      
#include "gpatInfo.hpp"

#if defined HAS_OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace vcflib;

struct copyNcounts{
  std::string       seqid       ;
  std::string       type        ;
  std::string       targetV     ;
  std::string       backgroundV ;
  std::string       callers     ;

  std::stringstream results ;
  long int          pos     ;
  long int          end     ;
  std::vector<double> target, background, total;
};

int nper = 1000;

void printHelp(void){
  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << "     vFst calculates vst, a measure of CNV stratification." << endl << endl;
 
  cerr << "Output : 4 columns :     "    << endl;
  cerr << "     1. seqid            "    << endl;
  cerr << "     2. position         "    << endl;
  cerr << "     3. end              "    << endl;
  cerr << "     3. vst              "    << endl;
  cerr << "     4. probability      "    << endl  << endl;

  cerr << "INFO: usage:  pVst --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf --type CN" << endl;
  cerr << endl;
  cerr << "INFO: required: t,target     -- argument: a zero based comma separated list of target individuals corresponding to VCF columns       "  << endl;
  cerr << "INFO: required: b,background -- argument: a zero based comma separated list of background individuals corresponding to VCF columns   "  << endl;
  cerr << "INFO: required: f,file       -- argument: a properly formatted VCF.                                                                  "  << endl;
  cerr << "INFO: required: y,type       -- argument: the genotype field with the copy number: e.g. CN|CNF                           "  << endl;
  cerr << "INFO: optional: r,region     -- argument: a tabix compliant genomic range : seqid or seqid:start-end                                 "  << endl;
  cerr << "INFO: optional: x,cpu        -- argument: number of CPUs [1] " << endl ;
  cerr << "INFO: optional: n,per        -- argument: number of permutations [1000] " << endl;       
  
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
//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of doubles

 Function does   : calculates the var

 Function returns: double

*/

inline double var(vector<double> & d, double mu){
  double variance = 0;

  for(vector<double>::iterator it = d.begin(); it != d.end(); it++){
    variance += pow((*it) - mu,2);
  }

  return variance / (d.size() - 1);
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of ints

 Function does   : calculates the mean

 Function returns: double

*/


inline double mean(vector<int> & data){

  double sum = 0;

  for(vector<int>::iterator it = data.begin(); it != data.end(); it++){
    sum += (*it);
  }
  return sum / data.size();
}

inline double mean(vector<double> & data){

  double sum = 0;

  for(vector<double>::iterator it = data.begin(); it != data.end(); it++){
    sum += (*it);
  }
  return sum / data.size();
}

double vst(copyNcounts * d){
  
  double muA = mean(d->total     );
  double vA = var(d->total, muA     );
  
  if(vA <= 0){
    return 0;
  }

  double muT = mean(d->target    );
  double muB = mean(d->background);

  double vT = var(d->target, muT    );
  double vB = var(d->background, muB);

  double cT = double(d->target.size())     / double(d->total.size());
  double cB = double(d->background.size()) / double(d->total.size());
  
  double ans =  (vA - ( (cT * vT) + (cB * vB) )) / vA;

  //purge negative
  if(ans < 0){
    ans = 0;
  }

  return ans;

}

void calc(copyNcounts * d){

  double v = vst(d);

  int nsuc   = 0;
  int trials = 0;

  for(uint i = 0; i < nper; i++){
    if (nsuc > 1){
      break;
    }
    
    trials += 1;

    std::random_shuffle(d->total.begin(), d->total.end());
    
    int tsize = d->target.size();
          
    d->target.clear();
    d->background.clear();
    
    int counter = 0;

    for(std::vector<double>::iterator it = d->total.begin();
	it!= d->total.end(); it++){
      if(counter < tsize){
	d->target.push_back(*it);
      }
      else{
	d->background.push_back(*it);
      }
      counter+=1;
    }
    
    //    std::cerr << "PER\t" << v << "\t" << vst(d) << std::endl;
    
    if (vst(d) >= v ){
      nsuc += 1;
    }
  }

  double p = 0;
  if(nsuc > 0){
    p = double(nsuc) / double(trials);
  }

  

  d->results << d->seqid 
	     << "\t" 
	     << d->pos 
	     << "\t" 
	     << d->end
	     << "\t"
	     << d->type
	     << "\t"
	     << v
	     << "\t"
	     << p ;
  
  


}



//------------------------------- SUBROUTINE --------------------------------

void loadIndices(map<int, int> & index, string set){
  
  vector<string>  indviduals = split(set, ",");

  vector<string>::iterator it = indviduals.begin();
  
  for(; it != indviduals.end(); it++){
    index[ atoi( (*it).c_str() ) ] = 1;
  }
}

// gotta load the dat so that we can permute using open MP
void loadDat(copyNcounts * d, 
	     std::string & type,
	     vector < map< string, vector<string> > > & target,
	     vector < map< string, vector<string> > > & background){

  for(vector < map< string, vector<string> > >::iterator it 
	= target.begin(); it!= target.end(); it++){
    d->target.push_back( atof((*it)[type].front().c_str()) );
    d->total.push_back(  atof((*it)[type].front().c_str()) );
    
    d->targetV += (*it)[type].front();
    d->targetV += ",";
    
  }
  for(vector < map< string, vector<string> > >::iterator it 
	= background.begin(); it!= background.end(); it++){
    d->background.push_back( atof((*it)[type].front().c_str()) );
    d->total.push_back( atof((*it)[type].front().c_str()) );  
  
    d->backgroundV += (*it)[type].front();
    d->backgroundV += ",";

  }

}

int main(int argc, char** argv) {

  std::srand ( unsigned ( std::time(0) ) );

  int cpu = 1;

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
	{"cpu"       , 1, 0, 'x'},
	{"per"       , 1, 0, 'n'},
	{"file"      , 1, 0, 'f'},
	{"target"    , 1, 0, 't'},
	{"background", 1, 0, 'b'},
	{"region"    , 1, 0, 'r'},
	{"type"      , 1, 0, 'y'},
  
	{0,0,0,0}
      };

    int index;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "n:r:d:t:b:f:y:x:hv", longopts, &index);
	switch (iarg)
	  {
	  case 'h':
	    {
	      printHelp();
	      return 0;
	    }
	  case 'n':
	    {
	      nper = atoi(optarg);
	      cerr << "INFO: pVst will do " << nper << "permutations." << endl;
	      break;
	    }
	  case 'v':
	    {
	      printVersion();
	      return 0;
	    }
	  case 'y':
	    {	    
	      type = optarg;
	      cerr << "INFO: Copy number will be found in : " << type << endl;
	      break;
	    }
	  case 'x':
	    {
	      cerr << "INFO: OpenMp will use " << optarg << " cpus. " << endl;
	      cpu = atoi(optarg);
	      break;
	    }
	  case 't':
	    {
	      loadIndices(ib, optarg);
	      cerr << "INFO: There are " << ib.size() << " individuals in the target" << endl;
	      cerr << "INFO: target ids: " << optarg << endl;
	      break;
	    }
	  case 'b':
	    {
	      loadIndices(it, optarg);
	      cerr << "INFO: There are " << it.size() << " individuals in the background" << endl;
	      cerr << "INFO: background ids: " << optarg << endl;
	      break;
	    }
	  case 'f':
	    {
	      cerr << "INFO: File: " << optarg  <<  endl;
	      filename = optarg;
	      break;
	    }
	  case 'd':
	    {
	      cerr << "INFO: only scoring sites where the allele frequency difference is greater than: " << optarg << endl;
	      deltaaf = optarg;
	      daf = atof(deltaaf.c_str());	    
	      break;
	    }
	  case 'r':
            cerr << "INFO: set seqid region to : " << optarg << endl;
	    region = optarg; 
	    break;
	  default:
	    break;
	  }
      }

#if defined HAS_OPENMP
    omp_set_num_threads(cpu);
#endif

    
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
 
    if(type == "NA"){
      cerr << "FATAL: failed to specify copy number genotype field" << endl;
      printHelp();
      return 1;
    }
 
    Variant var(variantFile);

    vector<string> samples = variantFile.sampleNames;
    int nsamples = samples.size();

    std::vector<copyNcounts*> dataBin;

    while (variantFile.getNextVariant(var)) {

      if(var.alt.size() > 1){
	continue;
      }

      std::map<string, bool> formatMap;
      for(std::vector<std::string>::iterator itz = var.format.begin();
	  itz != var.format.end(); itz++){
	formatMap[*itz] = true;
      }

      if(formatMap.find(type) == formatMap.end()){
	continue;
      }
      
 
      copyNcounts * varDat = new copyNcounts;

      varDat->pos   = var.position    ;
      varDat->seqid = var.sequenceName;
      varDat->type  = var.alt.front() ;
      
      if(var.info.find("CALLERS") != var.info.end()){
	stringstream caller;
	for(std::vector<std::string>::iterator z = var.info["CALLERS"].begin();
	    z != var.info["CALLERS"].end(); z++){
	  caller << (*z);
	  caller << ',';
	}
	
	varDat->callers = caller.str();
      }

      if(var.info.find("END") != var.info.end()){
	varDat->end = atol(var.info["END"].front().c_str());
      }
      else{
	varDat->end = varDat->pos;
      }

      vector < map< string, vector<string> > > target, background, total;
      
      int index = 0;
      
     
      for(int nsamp = 0; nsamp < nsamples; nsamp++){
	
	map<string, vector<string> > sample = var.samples[ samples[nsamp]];
	
	if(sample[type].front() != "."){
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
      
      loadDat(varDat, type, target, background);

      if(varDat->target.size() < 2 
	 || varDat->background.size() < 2){
	delete varDat;
	continue;
      }

      dataBin.push_back(varDat);

      // this odd pattern is for open MP ... later
      if(dataBin.size() == cpu){
	
#if defined HAS_OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
	for(int i = 0 ; i < dataBin.size(); i++){
	  calc(dataBin[i]);
	}
	for(int i = 0 ; i < dataBin.size(); i++){
	  std::cout << dataBin[i]->results.str() << "\t" << dataBin[i]->targetV << "\t" << dataBin[i]->backgroundV ;
	  
	  if(!dataBin[i]->callers.empty()){
	    std::cout << "\t" << dataBin[i]->callers << std::endl;
	  }
	  else{
	    std::cout << std::endl;
	  }
	  

	  delete dataBin[i];
	}
	dataBin.clear();
      }
    }
    return 0;		    
}
