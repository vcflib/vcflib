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
using namespace vcf;

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
  cerr << "      Splits VCF into unique files.  This resolves the issue " << endl;
  cerr << "      of same start. " << endl;

  cerr << "Output : N vcf files                 "    << endl;

  cerr << "INFO: usage:  splitUniqStarts --file my.vcf --prefix split " << endl;
  cerr << endl;
  cerr << "INFO: required: f,file    -- proper formatted VCF" << endl;
  cerr << "INFO: required: p,prefix    -- vcf prefix" << endl;


  printVersion();
}


int main(int argc, char** argv) {

  int maxNfiles = 1;

  string prefix;

  string filename;

  // set region to scaffold

  string region = "NA"; 

  // using vcflib; thanks to Erik Garrison 


  // zero based index for the target and background indivudals 
  

    const struct option longopts[] = 
      {
	{"version"   , 0, 0, 'v'},
	{"help"      , 0, 0, 'h'},
        {"file"      , 1, 0, 'f'},
	{"prefix"    , 1, 0, 'p'},
	{0,0,0,0}
      };

    int index;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "p:f:r:vh", longopts, &index);
	
	switch (iarg)
	  {
	  case 'p':
	    {
	      prefix = optarg;
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
	  default:
	    break;
	  }
	
      }
    
    if(filename.empty()){
      cerr << "FATAL: failed to specify a file" << endl;
      printHelp();
    }
    
    // scope since variant file has no close
    {
      VariantCallFile variantFile;
      
      if(!variantFile.open(filename)){
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
      
      Variant var(variantFile);
      
      std::map<long int, int> seen;
      
      long int lastP    = -1;
      std::string lastS = "NA";
     
      while (variantFile.getNextVariant(var)) {
	
	if(var.position < lastP && var.sequenceName == lastS){
	  
	  std::cerr << "FATAL: " << "VCF must be sorted." << std::endl;
	  std::cerr << "       lastPos = " << lastP << " currentPos = " << var.position << std::endl;
	  exit(1);
	}
	if(lastS != var.sequenceName){
	  for(std::map<long int, int>::iterator it = seen.begin();
	      it != seen.end(); it++){    
	    if(it->second > maxNfiles){
std::cerr << "INFO: " << lastS << "\t" << it->first << "\t" << it->second << std::endl;
	      maxNfiles = it->second;
	    }
	  }	
	  seen.clear();
	  lastS = var.sequenceName;
	}
	
	// counting non uniq starts
	if(var.position == lastP){
	  if(seen.find(var.position) != seen.end()){
	    seen[var.position] += 1;
	  }
	  else{
	    // start with 2 (current and last)
	    seen[var.position] = 2;
	  }
	}
	lastP = var.position;
      }
    }

    std::cerr << "INFO: " << maxNfiles << " required to uniq start." << std::endl;

    if(maxNfiles == 1){
      std::cerr << "INFO: all start positions are unique, bye!";
      return 0;
	
    }

    {
      VariantCallFile variantFile;
      if(!variantFile.open(filename)){
	std::cerr << "FATAL: could not open VCF" << endl;
	exit(1);
      }

      std::map<int, ofstream *> outputFiles;
      for(int i = 1; i <= maxNfiles; i++){
	ofstream * oz = new ofstream;
	
	stringstream fname;
        fname << prefix << "." << i << ".vcf";
	oz->open(fname.str().c_str());
	
	*oz << variantFile.header << std::endl;

	outputFiles[i] = oz;
      }

      Variant var(variantFile);
      
      int index = 1;

      while (variantFile.getNextVariant(var)) {
	*(outputFiles[index]) <<  var << std::endl;
	index += 1;
	if(index > maxNfiles){
	  index = 1;
	}
      }
      for(int i = 1; i <= maxNfiles; i++){
	outputFiles[i]->close();
      }
    }
    return 0;		    
}
