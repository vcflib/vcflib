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
  cerr << "     xpEHH estimates haplotype decay between the target and background populations.  Haplotypes are integrated                        " << endl;
  cerr << "     until EHH in the target and background is less than 0.05. The score is the itegrated EHH (target) / integrated EHH (background). " << endl;
  cerr << "     xpEHH does NOT integrate over genetic distance, as genetic maps are not availible for most non-model organisms.                  " << endl;

  cerr << "Output : 4 columns :      " << endl;
  cerr << "     1. seqid             " << endl;
  cerr << "     2. position          " << endl;
  cerr << "     3. alllele frequency " << endl;
  cerr << "     4. EHH-alt           " << endl;
  cerr << "     5. EHH-ref           " << endl;
  cerr << "     6. xpEHH             " << endl << endl;

  cerr << "INFO: xpEHH  --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf                                                " << endl;
  cerr << endl;

  cerr << "INFO: required: t,target     -- argument: a zero based comma separated list of target individuals corrisponding to VCF columns        " << endl;
  cerr << "INFO: required: b,background -- argument: a zero based comma separated list of background individuals corrisponding to VCF columns    " << endl;
  cerr << "INFO: required: f,file       -- argument: a properly formatted phased VCF file                                                        " << endl;
  cerr << "INFO: required: y,type       -- argument: type of genotype likelihood: PL, GL or GP                                                   " << endl;
  cerr << "INFO: optional: r,region     -- argument: a tabix compliant genomic range: seqid or seqid:start-end                                   " << endl;
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

void calc(string haplotypes[][2], int nhaps, vector<long int> pos, vector<double> afs, vector<int> & target, vector<int> & background, string seqid){

  for(int snp = 0; snp < haplotypes[0][0].length(); snp++){
    
    double ehhsat = 1;
    double ehhsab = 1;

    double ehhAT = 1 ;
    double ehhAB = 1 ;

    int start = snp;
    int end   = snp;
    int core  = snp; 

    while( ehhAT > 0.001 && ehhAB > 0.001 ) {
     
      start -= 1;
      end   += 1;
      
      if(start == -1){
	break;
      }
      if(end == haplotypes[0][0].length() - 1){
	break;
      }
      
      map<string , int> targetH, backgroundH;

      double sumrT = 0;
      double sumaT = 0;
      double sumrB = 0;
      double sumaB = 0;
      double nrefT = 0;
      double naltT = 0;
      double nrefB = 0;
      double naltB = 0;

      // hashing haplotypes into maps for both chr1 and chr2 target[][1 & 2]

      for(int i = 0; i < target.size(); i++){
	targetH[ haplotypes[target[i]][0].substr(start, (end - start)) ]++;
	targetH[ haplotypes[target[i]][1].substr(start, (end - start)) ]++;
      }
      for(int i = 0; i < background.size(); i++){
	backgroundH[ haplotypes[background[i]][0].substr(start, (end - start)) ]++;
	backgroundH[ haplotypes[background[i]][0].substr(start, (end - start)) ]++;
      }

      // interating over the target populations haplotypes
      for( map<string, int>::iterator th = targetH.begin(); th != targetH.end(); th++){    	
	// grabbing the core SNP (*th).first.substr((end-start)/2, 1) == "1"
	if( (*th).first.substr((end-start)/2, 1) == "1"){     
	   sumaT += r8_choose(th->second, 2);  
	   naltT += th->second;
	}
	else{
	  sumrT += r8_choose(th->second, 2);  
	  nrefT += th->second;
	}
      }

      // integrating over the background populations haplotypes
      for( map<string, int>::iterator bh = backgroundH.begin(); bh != backgroundH.end(); bh++){
	// grabbing the core SNP (*bh).first.substr((end-start)/2, 1) == "1"
        if( (*bh).first.substr((end - start)/2 , 1) == "1"){
	  sumaB += r8_choose(bh->second, 2);
	  naltB += bh->second;
        }
        else{
          sumrB += r8_choose(bh->second, 2);
          nrefB += bh->second;
        }
      }
      
      ehhAT = sumaT / (r8_choose(naltT, 2));
      ehhAB = sumaB / (r8_choose(naltB, 2));
            
      ehhsat += ehhAT;
      ehhsab += ehhAB;
    } 
    if(std::isnan(ehhsat) || std::isnan(ehhsab)){
      continue;
    }
    cout << seqid << "\t" << pos[snp] << "\t" << afs[snp] << "\t" << ehhsat << "\t" << ehhsab << "\t" << log(ehhsat/ehhsab) << endl;
  }   
}

double EHH(string haplotypes[][2], int nhaps){

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

  int phased = 0;

  string type = "NA";

    const struct option longopts[] = 
      {
	{"version"   , 0, 0, 'v'},
	{"help"      , 0, 0, 'h'},
        {"file"      , 1, 0, 'f'},
	{"target"    , 1, 0, 't'},
	{"background", 1, 0, 'b'},
	{"region"    , 1, 0, 'r'},
	{"type"      , 1, 0, 'y'},

	{0,0,0,0}
      };

    int findex;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "y:r:t:b:f:hv", longopts, &findex);
	
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
	  case 'b':
	    loadIndices(ib, optarg);
	    cerr << "INFO: there are " << ib.size() << " individuals in the background" << endl;
	    cerr << "INFO: background ids: " << optarg << endl;
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

      string sampleName = (*samp);
     
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
    vector<double>   afs;

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
	  calc(haplotypes, nsamples, positions, afs, target_h, background_h, currentSeqid);
	}
	clearHaplotypes(haplotypes, nsamples);
	positions.clear();
	currentSeqid = var.sequenceName;
	afs.clear();
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
      if(type == "GT"){
	populationTarget     = new gt();
        populationBackground = new gt();
        populationTotal      = new gt();
      }
      
      populationTarget->loadPop(target,         var.sequenceName, var.position);
      
      populationBackground->loadPop(background, var.sequenceName, var.position);
	
      populationTotal->loadPop(total,           var.sequenceName, var.position);
      
      
//      if(populationTotal->af > 0.99 || populationTotal->af < 0.01){
//	
//	delete populationTarget;
//	delete populationBackground;
//	delete populationTotal;
//	continue;
//      }

      afs.push_back(populationTotal->af);
      positions.push_back(var.position);
      loadPhased(haplotypes, populationTotal, nsamples);
      
      delete populationTarget;
      delete populationBackground;
      delete populationTotal;

    }

    calc(haplotypes, nsamples, positions, afs, target_h, background_h, currentSeqid);
    
    return 0;		    
}
