#include <iostream>
#include <fstream>
#include <getopt.h>
#include <map>
#include <list>
#include <vector>
#include <string>
#include "split.h"
#include <stdio.h> 
#include <stdlib.h>
#include "gpatInfo.hpp"
#include <math.h> 

using namespace std;

struct opts{
  bool     truncate;
  string   format;
  long int step;
  long int size;
  int      seqid;
  int      pos  ; 
  int      value;
};

struct score{
  long int position;
  double score;
};



void printHelp(void){
  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << "      Smoother averages a set of scores over a sliding genomic window.            " << endl;
  cerr << "      Smoother slides over genomic positions not the SNP indices. In other words  " << endl;
  cerr << "      the number of scores within a window will not be constant. The last         " << endl;
  cerr << "      window for each seqid can be smaller than the defined window size.          " << endl;
  cerr << "      Smoother automatically analyses different seqids separately.                " << endl;

  cerr << "Output : 4 columns :     "    << endl;
  cerr << "     1. seqid            "    << endl;
  cerr << "     2. window start     "    << endl;
  cerr << "     2. window end       "    << endl;
  cerr << "     3. averaged score   "    << endl  << endl;

  cerr << "INFO: usage: smoother --format pFst --file GPA.output.txt" << endl;
  cerr << endl;
  cerr << "INFO: required: f,file     -- argument: a file created by GPAT++                           " << endl;
  cerr << "INFO: required: o,format   -- argument: format of input file, case sensitive               " << endl;
  cerr << "                              available format options:                                    " << endl;
  cerr << "                                wcFst, pFst, bFst, iHS, xpEHH, abba-baba, col3             " << endl;
  cerr << "INFO: optional: w,window   -- argument: size of genomic window in base pairs (default 5000)" << endl;
  cerr << "INFO: optional: s,step     -- argument: window step size in base pairs (default 1000)      " << endl;
  cerr << "INFO: optional: t,truncate -- flag    : end last window at last position (zero based)      " << endl;
  printVersion();
  cerr << endl << endl;
}

double ngreater(list<score> & rangeData, double val){

  double n = 0;
 

  for(list<score>::iterator it = rangeData.begin(); 
      it != rangeData.end(); it++ ){
    if(it->score >= val){
      n += 1;
    }   
  }
  return n;
}


double windowAvg(list<score> & rangeData){

  double n = 0;
  double s = 0;

  for(list<score>::iterator it = rangeData.begin(); it != rangeData.end(); it++){
    s += it->score;
    n += 1;
  }
  return (s/n);
}

//calculation of Patterson's D statistic
double dStatistic(list<score> & rangeData){

  double abba = 0;
  double baba = 0;
  double dstat ;
  for(list<score>::iterator it = rangeData.begin(); it != rangeData.end(); it++){
    if(it->score == 0){ // means we have BABA locus
      baba += 1;
    }
    else{ // count towards ABBA locus
      abba += 1;
    }
  }
  dstat = (abba - baba) / (abba + baba ); // d-statistic implementation
  return (dstat);
}

void processSeqid(ifstream & file, string seqid, streampos offset, opts & opt){
    
  string line ;
  
  long int windowSize = opt.size;
  long int start = 0;
  long int end   = windowSize;

  list<score> windowDat;
  
  file.clear();
    
  file.seekg(offset);
  
  vector<string> sline;

  while(getline(file, line)){

    sline = split(line, '\t');     
    score current ;
    if(seqid != sline[opt.seqid]){
      break;
    }
    current.position = atol( sline[opt.pos].c_str() );
    current.score    = atof( sline[opt.value].c_str() );

    if(opt.format == "iHS"){
      current.score = fabs(current.score);
    }



    // add in if abba-baba to process second score. 


    if(current.position > end){

      double reportValue ;

      if(opt.format == "abba-baba"){
	reportValue = dStatistic(windowDat);
      }
      else{
	reportValue = windowAvg(windowDat);
      }
      // nan
      if( reportValue == reportValue){
	cout << seqid << "\t" << start << "\t" << end << "\t" << windowDat.size() << "\t" << reportValue;
	if(opt.format == "iHS"){
	  std::cout << "\t" << ngreater(windowDat, 2.5) ;
	}
	std::cout << std::endl;
      }
      
           
    }
    while(end < current.position){
      start += opt.step;
      end   += opt.step;
      while(windowDat.front().position < start && !windowDat.empty()){
	windowDat.pop_front();
      }
    }
    windowDat.push_back(current);  
  }
  // add function for D-stat if abba-baba
  double finalMean = windowAvg(windowDat);
  
  if(opt.truncate && (finalMean == finalMean) ){
    cout << seqid << "\t" << start << "\t" << windowDat.back().position - 1 << "\t" << windowDat.size()  << "\t" << finalMean;
  
    if(opt.format == "iHS"){
      std::cout << "\t" << ngreater(windowDat, 2.5) ;
    }
    std::cout << std::endl;

  }
  else if(finalMean == finalMean){
    cout << seqid << "\t" << start << "\t" << end << "\t" << windowDat.size() << "\t" << finalMean;

    if(opt.format == "iHS"){
      std::cout << "\t" << ngreater(windowDat, 2.5) ;
    }
    std::cout << std::endl;

  }
  cerr << "INFO: smoother finished : " << seqid << endl;
}

int main(int argc, char** argv) {
  
  map<string, int> acceptableFormats;
  acceptableFormats["pFst"]  = 1;
  acceptableFormats["col3"]  = 1;
  acceptableFormats["bFst"]  = 1;
  acceptableFormats["wcFst"] = 1;
  acceptableFormats["xpEHH"] = 1;
  acceptableFormats["iHS"]   = 1;
  acceptableFormats["cqf"]   = 1;
  acceptableFormats["deltaAf"]   = 1;
  acceptableFormats["abba-baba"]   = 1;

  opts opt;
  opt.size     = 5000;
  opt.step     = 1000;
  opt.format   = "NA";
  opt.truncate = false;

  string filename = "NA";

  static struct option longopts[] = 
    {
      {"version"   , 0, 0, 'v'},
      {"help"      , 0, 0, 'h'},
      {"file"      , 1, 0, 'f'},
      {"window"    , 1, 0, 'w'},
      {"step"      , 1, 0, 's'},
      {"format"    , 1, 0, 'o'},
      {"truncate"  , 0, 0, 't'},
      {0,0,0,0}
    };

  int index;
  int iarg=0;

  while(iarg != -1){
    iarg = getopt_long(argc, argv, "f:w:s:o:vht", longopts, &index);
    switch(iarg){
    case 't':
      {
	opt.truncate = true;
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
	filename = optarg;
	cerr << "INFO: file : " << filename << endl;
	break;
      }
    case 's':
      {
	opt.step = atol(optarg);
	cerr << "INFO: step size : " << optarg << endl;
	break;
      }
    case 'w':
      {
	opt.size = atol(optarg);
	cerr << "INFO: step size : " << optarg << endl;
	break;
      }
    case 'o':
      {
	opt.format = optarg;
	cerr << "INFO: specified input format : " << optarg << endl;
	break;
      }
    }
  }
  if(filename == "NA"){
    cerr << "FATAL: file was not specified!" << endl << endl;
    printHelp();
    return 1;
  }

  if(acceptableFormats.find(opt.format) == acceptableFormats.end()){
    cerr << "FATAL: unacceptable input file format, see --format "  << endl << endl;
    printHelp();
    return 1;
  }
  if(opt.format == "deltaAf"){
    opt.seqid = 0;
    opt.pos   = 1;
    opt.value = 4;
  }
  else if(opt.format == "abba-baba"){
    opt.seqid = 0;
    opt.pos   = 1;
    opt.value = 2;
  }
  else if(opt.format == "col3"){
    opt.seqid = 0;
    opt.pos   = 1;
    opt.value = 2;
  }
  else if(opt.format == "pFst"){
    opt.seqid = 0;
    opt.pos   = 1;
    opt.value = 2;
  }
  else if (opt.format == "bFst"){
    opt.seqid = 0;
    opt.pos   = 1;
    opt.value = 8;
  }
  else if (opt.format == "wcFst"){
    opt.seqid = 0;
    opt.pos   = 1;
    opt.value = 4;
  }
  else if(opt.format == "cqf"){
    opt.seqid = 0;
    opt.pos   = 1;
    opt.value = 5;
  }
  else if(opt.format == "xpEHH"){
    opt.seqid = 0;
    opt.pos   = 1;
    opt.value = 5;
  }
  else if(opt.format == "iHS"){
    opt.seqid = 0;
    opt.pos   = 1;
    opt.value = 6;
  }

  else{
    cerr << "FATAL: input format flag not specified correctly : " << opt.format << endl;
    cerr << "INFO:  please use smoother --help" << endl;
    return 1;
  }
  
  ifstream ifs(filename.c_str());
 
  string currentSeqid = "NA";

  string line;

  map<string, streampos > seqidIndex;
  
  if(ifs){
    while(getline(ifs, line)){
      vector<string> sline = split(line, '\t');
      if(sline[opt.seqid] != currentSeqid){
	
	long int bline = ifs.tellg() ;
	bline -=  ( line.size() +1 );
	
	//	std::cerr << "INFO: seqid: " << sline[opt.seqid] << " tellg: " << bline << std::endl;

	map<string, streampos>::iterator it;

	if(seqidIndex.find(sline[opt.seqid]) != seqidIndex.end() ){
	  cerr << "FATAL: file is unsorted!" << endl;
	  return 1;
	}
	seqidIndex[sline[opt.seqid]] = bline;
	currentSeqid = sline[opt.seqid];
      }
    }
  }
  else{
    cerr << "FATAL: no lines -- or -- couldn't open file" << endl;
  }

  for( map< string, streampos>::iterator it = seqidIndex.begin(); it != seqidIndex.end(); it++){
    cerr << "INFO: processing seqid : "<< (it->first) << endl;
    processSeqid(ifs, (it->first),(it->second), opt);
  }
  
  ifs.close();
  cerr << "INFO: smoother has successfully finished" << endl;

  return 0;

}
