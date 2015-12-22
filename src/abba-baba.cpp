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
  cerr << "     abba-baba calculates the tree pattern for four indviduals.                         " << endl;
  cerr << "     This tool assumes reference is ancestral and ignores non abba-baba sites.          " << endl;
  cerr << "     The output is a boolian value: 1 = true , 0 = false for abba and baba.             " << endl;
  cerr << "     the tree argument should be specified from the most basal taxa to the most derived." << endl;
  cerr <<  endl;
  cerr << "     Example:" << endl;

  cerr <<  "     D   C  B   A "  << endl;
  cerr <<  "     \\ /  /    /  "  << endl;
  cerr <<  "      \\  /    /   "  << endl;
  cerr <<  "       \\    /    "  << endl;
  cerr <<  "        \\  /     "  << endl;
  cerr <<  "         /        "  << endl;
  cerr <<  "        /         "  << endl;
  
  cerr << " --tree A,B,C,D"  << endl << endl;
    

  cerr << "Output : 4 columns :     "    << endl;
  cerr << "     1. seqid            "    << endl;
  cerr << "     2. position         "    << endl;
  cerr << "     3. abba             "    << endl;
  cerr << "     4. baba             "    << endl;

  cerr << "INFO: usage:  abba-baba --tree 0,1,2,3 --file my.vcf --type PL" << endl;
  cerr << endl;
  cerr << "INFO: required: t,tree       -- a zero based comma seperated list of target individuals corrisponding to VCF columns" << endl;
  cerr << "INFO: required: f,file       -- a properly formatted VCF.                                                           " << endl;
  cerr << "INFO: required: y,type       -- genotype likelihood format ; genotypes: GP,GL or PL;                                " << endl;
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

/*random sample heterozygous genotypes could eventually be weighted 
by genotype likelihoods  and added complexity for linked, phased genos
random sampling adds noise but will not affect the overall measurement
of D-statistic */
int  sample_het(int &rv){
  rv = rand() % 2 ; // pick from 0/1 het with 50-50 odds
}


int  containsAlt(string gt){
  if(gt == "1/1"){
    return 1;
  }
  if(gt == "1|1"){
    return 1;
  }
  // heterozygous cases need to be randomly sampled for diploids
  int rv = 0 ;
  if(gt == "0/1"){
    return rv = sample_het(rv);
  }
  if(gt == "0|1"){
    return rv = sample_het(rv);
  }
  if(gt == "1|0"){
    return rv = sample_het(rv);
  }
  // all else return zero state
  return 0;
}


// flag as unused func'n
int  containsRef(string gt){
  if(gt == "0/1"){
    return 1;
  }
  if(gt == "0/0"){
    return 1;
  }
  if(gt == "0|0"){
    return 1;
  }

  if(gt == "0|1"){
    return 1;
  }
  if(gt == "1|0"){
    return 1;
  }
  return 0;
}

void loadIndices(vector<int> & tree, string set){
  
  vector<string>  indviduals = split(set, ",");

  if(indviduals.size() < 4){
    cerr << "FATAL: the abba-baba requires four indviduals provided to the tree option" << endl;
    exit(1);
  }
  for( vector<string>::iterator it = indviduals.begin(); it != indviduals.end(); it++){
    
    int indx = atoi((*it).c_str());
    cerr << indx << endl; //print sample index for user check
    tree.push_back(indx);
  }
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

  // zero based index for the tree
  
  vector<int> tree;
  
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
	{"tree"      , 1, 0, 't'},
	{"background", 1, 0, 'b'},
	{"region"    , 1, 0, 'r'},
	{"type"      , 1, 0, 'y'},
	{0,0,0,0}
      };

    int index;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "r:d:t:f:y:hv", longopts, &index);
	
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
	    break;
	  case 't':
	    loadIndices(tree, optarg);
	    break;
	  case 'f':
	    cerr << "INFO: File: " << optarg  <<  endl;
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

    
    if(filename == "NA"){
      cerr << "FATAL: did not specify the file\n";
      printHelp();
      exit(1);
    }
    
    variantFile.open(filename);
    
    if(region != "NA"){
      variantFile.setRegion(region);
    }

    if (!variantFile.is_open()) {
      exit(1);
    }
    map<string, int> okayGenotypeLikelihoods;

    okayGenotypeLikelihoods["PL"] = 1;
    okayGenotypeLikelihoods["PO"] = 1;
    okayGenotypeLikelihoods["GL"] = 1;
    okayGenotypeLikelihoods["GP"] = 1;


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

    Variant var(variantFile);
    
    vector<string> sampleNames = variantFile.sampleNames;

    srand(time(0)); //initialize random number generator

    while (variantFile.getNextVariant(var)) {

      if(var.alt.size() > 1){
	continue;
      }
      
      map<string, vector<string> >  tA, tB, tC, tD;

      tA = var.samples[sampleNames[tree[0]]];
      tB = var.samples[sampleNames[tree[1]]];
      tC = var.samples[sampleNames[tree[2]]];
      tD = var.samples[sampleNames[tree[3]]];

      if(tA["GT"].front() == "./." || tB["GT"].front() == "./." || tC["GT"].front() == "./." || tD["GT"].front() == "./."){
	continue;
      }
      
      int A = 0,B = 0,C = 0,D = 0; // set default allelic state to zero

      double abba = 0; //booleans for abab or baba state.
      double baba = 0;

      A = containsAlt(tA["GT"].front());
      B = containsAlt(tB["GT"].front());
      C = containsAlt(tC["GT"].front());
      D = containsAlt(tD["GT"].front());

      if(D == 1 && C == 0 && B == 0 && A == 1){
	abba = 1;
      }
      if(D == 0 && C == 1 && B == 0 && A == 1){
	baba = 1;
      }

      if(D == 0 && C == 1 && B == 1 && A == 0){
        abba = 1;
      }
      if(D == 1 && C == 0 && B == 1 && A == 0){
	baba = 1;
      }


      if(abba == 0 && baba == 0){
	continue;
      }

      cout << var.sequenceName << "\t" << var.position << "\t" << abba << "\t" << baba << endl;
      //cout << var.sequenceName << "\t" << var.position << "\t" << abba << "\t" << baba << "\t" << A << B << C << D << endl;
      // above is alternate print to check that we are getting observed
      // ABBA or BABA patterns
    }
    return 0;		    
}
