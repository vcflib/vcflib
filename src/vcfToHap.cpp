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
// maaas speed
#include <omp.h>
// print lock
omp_lock_t lock;



struct opts{
  int         threads             ;
  std::string filename            ;
  std::string mapFile             ;
  std::string seqid               ;
  std::string geneticMapFile      ;
  std::string type                ;
  std::string region              ;
  std::map<int, double> geneticMap;
  double      af                  ;

}globalOpts;


using namespace std;
using namespace vcflib;

void printHelp(void){
  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << "     iHS calculates the integrated ratio of haplotype decay between the reference and non-reference allele. " << endl;
  

  cerr << "Output : 4 columns :                  "    << endl;
  cerr << "     1. seqid                         "    << endl;
  cerr << "     2. position                      "    << endl;
  cerr << "     3. target allele frequency       "    << endl;
  cerr << "     4. integrated EHH (alternative)  "    << endl;
  cerr << "     5. integrated EHH (reference)    "    << endl;
  cerr << "     6. iHS ln(iEHHalt/iEHHref)      "    << endl  << endl;

  cerr << "Usage:" << endl;

  cerr << "      iHS  --target 0,1,2,3,4,5,6,7 --file my.phased.vcf  \\" << endl; 
  cerr << "           --region chr1:1-1000 > STDOUT 2> STDERR          " << endl << endl;

  cerr << "Params:" << endl;
  cerr << "       required: t,target  <STRING>  A zero base comma separated list of target" << endl;
  cerr << "                                     individuals corrisponding to VCF columns  " << endl;
  cerr << "       required: r,region  <STRING>  A tabix compliant genomic range           " << endl;
  cerr << "                                     format: \"seqid:start-end\" or \"seqid\"  " << endl; 
  cerr << "       required: f,file    <STRING>  Proper formatted and phased VCF.          " << endl;
  cerr << "       required: y,type    <STRING>  Genotype likelihood format: GT,PL,GL,GP   " << endl;
  cerr << "       optional: a,af      <DOUBLE>  Alterantive alleles with frquences less   " << endl; 
  cerr << "                                     than [0.05] are skipped.                  " << endl;
  cerr << "       optional: x,threads <INT>     Number of CPUS [1].                       " << endl;
  cerr << "       recommended: g,gen <STRING>   A PLINK formatted map file.               " << endl;
  cerr << endl;
 
  printVersion();

  exit(1);
}


bool gDist(int start, int end, double * gd){
  
  if(globalOpts.geneticMap.find(start) == globalOpts.geneticMap.end()){
    return false;
  }
  if(globalOpts.geneticMap.find(end) == globalOpts.geneticMap.end()){
    return false;
  }
  *gd = abs(globalOpts.geneticMap[start] - globalOpts.geneticMap[end]);
  return true;
}

void loadGeneticMap(int start, int end){

  if(globalOpts.geneticMapFile.empty()){
    std::cerr << "WARNING: No genetic map." << std::endl;
    std::cerr << "WARNING: A constant genetic distance is being used: 0.001." << std::endl;
    return;
  }

  ifstream featureFile (globalOpts.geneticMapFile.c_str());
 
  string line;

  int lastpos      = 0;
  double lastvalue = 0;

  if(featureFile.is_open()){

    while(getline(featureFile, line)){

      vector<string> region = split(line, "\t");

      if(region.front() != globalOpts.seqid){
	std::cerr << "WARNING: seqid MisMatch: " << region.front() << " " << globalOpts.seqid << std::endl;
	continue;
      }

      int   pos = atoi(region[3].c_str()) ;
      double cm = atof(region[2].c_str()) ;

      if(lastpos == 0 && start > pos){
	lastpos = pos;
	continue;
      }
     


      int diff     = abs(pos - lastpos);
      double vdiff = abs(lastvalue - cm );
      double chunk = vdiff/double(diff);

      //      std::cerr << "INFO: " << diff << " " << vdiff << endl;

      double running = lastvalue;

      for(int i = lastpos; i < pos; i++){
	globalOpts.geneticMap[i] = running;
	running += chunk;
      }

      if(pos > end){
	break;
      }


      lastpos = pos;
      lastvalue = cm;
    }
  }

  featureFile.close();

  if(globalOpts.geneticMap.size() < 1){
    std::cerr << "FATAL: Problem loading genetic map" << std::endl;
    exit(1);
  }
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

  globalOpts.threads = 1   ;
  globalOpts.af      = 0.05;

  // zero based index for the target and background indivudals 
  
  map<int, int> it, ib;
  
    const struct option longopts[] = 
      {
	{"version"   , 0, 0, 'v'},
	{"help"      , 0, 0, 'h'},
        {"file"      , 1, 0, 'f'},
	{"target"    , 1, 0, 't'},
	{"region"    , 1, 0, 'r'},
	{"gen"       , 1, 0, 'g'},
	{"type"      , 1, 0, 'y'},
	{"threads"   , 1, 0, 'x'},
	{"af"        , 1, 0, 'a'},
	{0,0,0,0}
      };

    int findex;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "a:x:g:y:r:d:t:b:f:hv", longopts, &findex);
	
	switch (iarg)
	  {
	  case 'a':
	    {
	      globalOpts.af = atof(optarg);
	      break;
	    }
	  case 'x':
	    {
	      globalOpts.threads = atoi(optarg);
	      break;
	    }
	  case 'g':
	    {
	      globalOpts.geneticMapFile = optarg;
	      break;
	    }
	  case 'h':
	    {
	      printHelp();
	      break;
	    }
	  case 'v':
	    {
	      printVersion();
	      break;
	    }
	  case 'y':
	    {
	      globalOpts.type = optarg;
	      break;
	    }
	  case 't':
	    {
	      loadIndices(it, optarg);
	      cerr << "INFO: there are " << it.size() << " individuals in the target" << endl;
	      cerr << "INFO: target ids: " << optarg << endl;
	      break;
	    }
	  case 'f':
	    {
	      cerr << "INFO: file: " << optarg  <<  endl;
	      globalOpts.filename = optarg;
	      break;
	    }
	  case 'r':
	    {
	      cerr << "INFO: set seqid region to : " << optarg << endl;
	      globalOpts.region = optarg; 
	      break;
	    default:
	      break;
	    }
	  }
      }

  omp_set_num_threads(globalOpts.threads);

    map<string, int> okayGenotypeLikelihoods;
    okayGenotypeLikelihoods["PL"] = 1;
    okayGenotypeLikelihoods["GL"] = 1;
    okayGenotypeLikelihoods["GP"] = 1;
    okayGenotypeLikelihoods["GT"] = 1;
    

    // add an option for dumping

//    for(std::map<int, double>::iterator gm = geneticMap.begin(); gm != geneticMap.end(); gm++){
//      cerr << "pos: " << gm->first << " cm: " << gm->second << endl; 
//    }

    if(globalOpts.type.empty()){
      cerr << "FATAL: failed to specify genotype likelihood format : PL or GL" << endl;
      printHelp();
      exit(1);
    }
    if(okayGenotypeLikelihoods.find(globalOpts.type) == okayGenotypeLikelihoods.end()){
      cerr << "FATAL: genotype likelihood is incorrectly formatted, only use: PL or GL" << endl;
      printHelp();
      exit(1);
    }

    if(globalOpts.filename.empty()){
      cerr << "FATAL: did not specify a file" << endl;
      printHelp();
      exit(1);
    }

    if(it.size() < 2){
      cerr << "FATAL: target option is required -- or -- less than two individuals in target\n";
      printHelp();
      exit(1);
    }

    // using vcflib; thanksErik 

    VariantCallFile variantFile;

    variantFile.open(globalOpts.filename);
    
    if(globalOpts.region.empty()){
      cerr << "FATAL: region required" << endl;
      exit(1);
    }
    if(! variantFile.setRegion(globalOpts.region)){
      cerr <<"FATAL: unable to set region" << endl;
      exit(1);
    }

    if (!variantFile.is_open()) {
      exit(1);
    }
    
    Variant var( variantFile );
    vector<int> target_h, background_h;

    int index   = 0; 
    int indexi  = 0;


    vector<string> samples = variantFile.sampleNames;
    int nsamples = samples.size();

    for(vector<string>::iterator samp = samples.begin(); samp != samples.end(); samp++){
      
      string sampleName = (*samp);
     
      if(it.find(index) != it.end() ){
	target_h.push_back(indexi);
	indexi++;
      }
      index++;
    }
    
   
    vector<long int> positions;
    
    vector<double> afs;

    string haplotypes [target_h.size()][2];    
    

    while (variantFile.getNextVariant(var)) {

      globalOpts.seqid = var.sequenceName;

      if(!var.isPhased()){
	cerr << "FATAL: Found an unphased variant. All genotypes must be phased!" << endl;
	exit(1);
      }

      if(var.alleles.size() > 2){
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
      
      if(globalOpts.type == "PL"){
	populationTarget     = new pl();
      }
      if(globalOpts.type == "GL"){
	populationTarget     = new gl();
      }
      if(globalOpts.type == "GP"){
	populationTarget     = new gp();
      }
      if(globalOpts.type == "GT"){
	populationTarget     = new gt();
      }

      populationTarget->loadPop(target, var.sequenceName, var.position);
      
      if(populationTarget->af <= globalOpts.af || populationTarget->af >= (1-globalOpts.af) ){
	delete populationTarget;
	continue;
      }
      positions.push_back(var.position);
      afs.push_back(populationTarget->af);
      loadPhased(haplotypes, populationTarget, populationTarget->gts.size()); 
    
      populationTarget = NULL;
      delete populationTarget;
    }

    if(!globalOpts.geneticMapFile.empty()){
      cerr << "INFO: loading genetics map" << endl;
      loadGeneticMap(positions.front(), positions.back());
      cerr << "INFO: finished loading genetics map" << endl;
    }
    
    if(positions.size() != haplotypes[0][0].size()){
      std::cerr << "FATAL: there are " << positions.size() << " and " << haplotypes[0][0].size() << " haplotype SNPs" << endl;
      exit(1);
    }


    for(int i = 0; i < positions.size(); i++){
      cerr.precision(8);

      std::cerr << "rs" << i << " " << positions[i] << " " << globalOpts.geneticMap[positions[i]] << " " << "A T" << endl;
    }
    
    for(int i = 0; i < target_h.size() ; i++){
      cout << haplotypes[i][0][0] ;
      for(int j = 1; j < haplotypes[i][1].size(); j++){
	cout << " " << haplotypes[i][0][j] ;
      }
      cout << endl;
      cout << haplotypes[i][1][0] ;
      for(int j = 1; j < haplotypes[i][1].size(); j++){
        cout << " " << haplotypes[i][1][j] ;
      }
      cout << endl;
      
    }



    clearHaplotypes(haplotypes, target_h.size());

    exit(0);		    

}
