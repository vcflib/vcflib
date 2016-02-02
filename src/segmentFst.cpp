
/*

This program was created at:  Tue Sep  8 21:05:23 2015
This program was created by:  Zev N. Kronenberg


Contact: zev.kronenber@gmail.com

Organization: Unviersity of Utah
    School of Medicine
    Salt Lake City, Utah


The MIT License (MIT)

Copyright (c) <2015> <Zev N. Kronenberg>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


*/

#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <algorithm>
#include "split.h"
#include "gpatInfo.hpp"

struct options{
  std::string file;
  double cut;
}globalOpts;

using namespace std;

static const char *optString = "hf:s:";

void printHelp(void){

  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << "      Creates genomic segments (bed file) for regions with high wcFst  " << endl;

  cerr << "Output : 8 columns :                 "    << endl;
  cerr << "     1. Seqid                        "    << endl;
  cerr << "     2. Start (zero based)           "    << endl;
  cerr << "     3. End   (zero based)           "    << endl;
  cerr << "     4. Average Fst                  "    << endl;
  cerr << "     5. Average high Fst (Fst > -s)  "    << endl;
  cerr << "     6. N Fst values in segment      "    << endl;
  cerr << "     7. N high fst values in segment "    << endl;
  cerr << "     8. Segment length               "    << endl;

  cerr << "INFO: usage:  segmentFst -s 0.7 -f wcFst.output.txt " << endl;
  cerr << endl;
  cerr << "INFO: required: -f            -- Output from wcFst     "   << endl;
  cerr << "INFO: optional: -s            -- High Fst cutoff [0.8] "    << endl;

  cerr << endl;

  printVersion();
}

//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
    {
    int opt = 0;
    globalOpts.file = "NA";
    opt = getopt(argc, argv, optString);
    while(opt != -1){
      switch(opt){
      case 's':
	{
	  string op = optarg;    
	  globalOpts.cut = atof(op.c_str());
	  break;
	}
      case 'h':
	{
	  printHelp();
	  exit(1);
	  break;
	}
	
      case 'f':
	{
	  globalOpts.file = optarg;
	  break;
	}
      case '?':
	{
	  break;
	}
      }	
	opt = getopt( argc, argv, optString ); 
   }
return 1;
}
//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  :

 Function does   :

 Function returns:

*/
void process(vector<int> & pos, vector<double> & value, vector<string> & seqid)
{
  if(value.size() < 10){
    return;
  }

  vector<int> yes;

  vector<double> sortedVs;
  
  for(int i = 0; i < value.size(); i++){
    yes.push_back(0);
    sortedVs.push_back(value[i]);
  }

  sort(sortedVs.begin(), sortedVs.end());
  
  int n = value.size()-1;

  while( (double(n) / double(value.size())) > 0.99){
    n-=1;
  }

  double high = globalOpts.cut;

  cerr << "95th qunatile: " << high << " " << endl;
  
  for(int i = 1; i < 10000; i = (i*2)){
    if(i >= value.size()){
      break;
    }
    for(int j= 0 ; j < value.size() -i ; j++){
      
      int nbelow = 0;
      int nabove = 0;
      
      for(int z = j; z < i+j; z++){
	if(value[z] > high){
	  nabove +=1;
	}
	else{
	  nbelow += 1;
	}
      }
      if((2*nabove) > nbelow){
	for(int fl = j; fl < i+j; fl++){
	  yes[fl] += 1;
	}
      }
    }
  }
  
  for(int y = 0; y < yes.size(); y++){
    int begin = y;
    int end   = y;
    if(yes[y] > 0){
      while(yes[y] > 0){
	end = y;
	y  += 1;
      }

      double fstSumHigh = 0;
      double nHigh  = 0;
      double fstSum = 0;
      double n      = 0;
      for(int t = begin; t <= end; t++){
	
	if(value[t] > high){
	  nHigh += 1;
	  fstSumHigh += value[t];
	}
	
	n += 1;
	fstSum += value[t];
      }

      double avgFst = (fstSum / n);
      double avgFstHigh = fstSumHigh / nHigh;

      cout << seqid[begin] << "\t" 
	   << pos[begin] -1 << "\t" 
	   << pos[end] - 1  << "\t" 
	   << avgFst << "\t"
	   << avgFstHigh << "\t"
	   << n << "\t"
	   << nHigh << "\t"
	   << (pos[end] - pos[begin]) 
	   << endl;
    }
  }
}

//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
  globalOpts.cut = 0.8;
  int parse = parseOpts(argc, argv);

  string last = "NA";

  vector<string> seqid;
  vector<int>      pos;
  vector<double> value;

  if(globalOpts.file.empty()){
    printHelp();
    exit(1);
  }

  string line;
  ifstream myfile (globalOpts.file);
  if (myfile.is_open())
    {
      while ( getline (myfile,line) )
	{
	  vector<string> lineDat = split(line, '\t');
	  if(lineDat.size() != 5){
	    cerr << "FATAL: not valid input" << endl;
	    exit(1);
	  }
	  if(last.compare(lineDat[0]) != 0){
	    last = lineDat[0];
	    process(pos, value, seqid);
	    pos.clear();
	    value.clear();
	    seqid.clear();
	  }
	  else{
	    // WARNING must have af greater than 0.05 <- real variant ....

	    if(atof(lineDat[2].c_str()) > 0.05 || atof(lineDat[3].c_str()) > 0.05){
	      seqid.push_back(lineDat[0]);
	      pos.push_back(atoi(lineDat[1].c_str()));
	      value.push_back(atof(lineDat[4].c_str()));	  
	    }
	  }
	}
      process(pos, value, seqid);
      myfile.close();
    }
  return 0;
}
