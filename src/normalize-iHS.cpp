
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
  double afDiff;
}globalOpts;

struct iHSdat{
  std::string seqid;
  std::string start;
  std::string F1   ;
  std::string F2   ;

  double af   ;
  double ehhR ;
  double ehhA ;
  double iHS  ;
  double niHS ;
};


using namespace std;

static const char *optString = "hf:s:";

void printHelp(void){

  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << "      normalizes iHS or XP-EHH scores  " << endl;

  cerr << "Output : normalize-iHS adds one additional column to input (normalized score)." << endl;

  cerr << "INFO: usage:  normalizeHS -s 0.01 -f input.txt " << endl;
  cerr << endl;
  cerr << "INFO: required: -f            -- Output from iHS or XPEHH "   << endl;
  cerr << "INFO: optional: -s            -- Max AF diff for window [0.01]"   << endl;

  cerr << endl;

  printVersion();
}

//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
    {
    int opt = 0;
    opt = getopt(argc, argv, optString);
    while(opt != -1){
      switch(opt){
      case 's':
	{
	  string op = optarg;    
	  globalOpts.afDiff = atof(op.c_str());
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

bool sortAF(iHSdat * L, iHSdat * R){
  if(L->af < R->af){
    return true;
  }
  return false;
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of doubles

 Function does   : calculates the var

 Function returns: double

*/

double var(vector<double> & data, double mu){
  double variance = 0;

  for(vector<double>::iterator it = data.begin(); it != data.end(); it++){
    variance += pow((*it) - mu,2);
  }

  return variance / (data.size() - 1);
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of doubles

 Function does   : computes the mean

 Function returns: the mean

*/


double windowAvg(std::vector<double> & rangeData){

  long double n = 0;
  long double s = 0;

  for(std::vector<double>::iterator it = rangeData.begin(); it != rangeData.end(); it++){
    s += *it;
    n += 1;
  }
  

  return (s/n);
}



//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of iHS data 

 Function does   : normalizes

 Function returns: nothing

*/

void normalize(std::vector<iHSdat *> & data, int * pos){
  
  std::vector<double> windat;

  int start = *pos;
  int end   = *pos;

  while((abs(data[start]->af - data[end]->af ) < globalOpts.afDiff) 
	&& end < data.size() -1 ){
    end += 1;
  }

  for(int i = start; i <= end; i++){
    windat.push_back(data[i]->iHS);
  }
  
  double avg = windowAvg(windat);
  double sd  = sqrt(var(windat, avg));

  std::cerr << "start: " << data[start]->af << " " 
	    << "end: " << data[end]->af  << " " 
            << "n iHS scores: " << windat.size() << " "  
	    << "mean: " << avg << " " 
	    << "sd: " << sd << std::endl;

  for(int i = start; i <= end; i++){
    data[i]->niHS = (data[i]->iHS - avg) / (sd);
  }
  
  *pos = end;
  
}

//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
  globalOpts.afDiff = 0.01;
  int parse = parseOpts(argc, argv);

  if(globalOpts.file.empty()){
    std::cerr << "FATAL: no file" << std::endl;
    exit(1);
  }

  std::vector<iHSdat *> data;

  string line;
  ifstream myfile (globalOpts.file);
  if (myfile.is_open())
    {
      while ( getline (myfile,line) ){
	vector<string> lineDat = split(line, '\t');

	iHSdat * tp = new iHSdat;
	tp->seqid = lineDat[0];
	tp->start = lineDat[1];
	tp->af    = atof(lineDat[2].c_str());
	tp->ehhR  = atof(lineDat[3].c_str());
	tp->ehhA  = atof(lineDat[4].c_str());
	tp->iHS   = atof(lineDat[5].c_str());
	tp->F1    = lineDat[6].c_str();
	tp->F2    = lineDat[7].c_str();
	tp->niHS  = 0;

	data.push_back(tp);

      }
  
    myfile.close();
    }
  else{
    cerr << "FATAL: could not open file: " << globalOpts.file << endl;
    exit(1);
  }
  

  std::cerr << "INFO: sorting " << data.size() << " scores by AF" << std::endl;

  sort(data.begin(), data.end(), sortAF);

  std::cerr << "INFO: finished sorting" << std::endl;

  for(int i = 0; i < data.size() ; i++){
    normalize(data, &i);
  }


  for(int i = 0; i < data.size(); i++){
    std::cout << data[i]->seqid << "\t"
	      << data[i]->start << "\t"
	      << data[i]->af << "\t"
	      << data[i]->ehhR << "\t"
	      << data[i]->ehhA << "\t"
	      << data[i]->iHS << "\t"
	      << data[i]->niHS << "\t"
	      << data[i]->F1 << "\t"
	      << data[i]->F2 << std::endl;
	     
  }


  return 0;
}
