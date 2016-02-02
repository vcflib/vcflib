/*

This program was created at:  Fri Apr 17 14:59:53 2015
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
#include <fstream>
#include "split.h"
#include <vector>
#include <string>
#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "gpatInfo.hpp"

struct options{
  std::string file;
  int npermutation;
  int nsuc; 
}globalOpts;

static const char *optString = "f:n:s:";

using namespace std;


//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
    {
    int opt = 0;
    globalOpts.file = "NA";

    globalOpts.nsuc         = 1;
    globalOpts.npermutation = 1000;
    
    opt = getopt(argc, argv, optString);
    while(opt != -1){
      switch(opt){
      case 'f':
	{
	  globalOpts.file =  optarg;
	  break;
	}
      case 'n':
	{
	  globalOpts.npermutation = atoi(((string)optarg).c_str());
	  cerr << "INFO: permuteGPAT++ will do N permutations: " << globalOpts.npermutation << endl;
	  break;
	}
      case 's':
	{
	  globalOpts.nsuc = atoi(((string)optarg).c_str());
	  cerr << "INFO: permuteGPAT++ will stop permutations after N successes: " << globalOpts.nsuc << endl;
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
 Function input  : NA

 Function does   : prints help

 Function returns: NA

*/
void printHelp()
{

  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << "     permuteGPAT++ is a method for adding empirical p-values to a GPAT++ score." << endl ;
  cerr << "     Currently permuteGPAT++ only supports wcFst, but will be extended.    " << endl ;
  cerr << endl;
  cerr << "OUTPUT: permuteGPAT++ will append three additional columns:" << endl;
  cerr << "        1. The number of successes                         " << endl;
  cerr << "        2. The number of trials                            " << endl;
  cerr << "        3. The empirical p-value                           " << endl << endl; 

  cerr << "INFO: usage:  permuteGPAT++ -f gpat.txt -n 5 -s 1 "<< endl;
  cerr << endl;
  cerr << "INFO: file:    f   -- argument: the input file     "<< endl;
  cerr << "INFO: number:  n   -- argument: the number of permutations to run for each value [1000]" << endl;
  cerr << "INFO: success: s   -- argument: stop permutations after \'s\' successes [1]"             << endl;

  cerr << endl;
  printVersion();

}

//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
int parse = parseOpts(argc, argv);

 if(globalOpts.file.compare("NA") == 0){
   cerr << "FATAL: no file was provided" << endl;
   printHelp();
   exit(1);
 }


 vector<double> data;

 ifstream gpat (globalOpts.file.c_str());

 string line;

 if(gpat.is_open()){

   while(getline(gpat, line)){
     vector<string> region = split(line, "\t");
     // will change for other output
     double fst = atof(region[4].c_str());

     if(fst < 0){
       fst = 0;
     }

     data.push_back(fst);
   }
 }
 else{
   cerr << "FATAL: coult not open file: " << globalOpts.file << endl;
   exit(1);
 }

 gpat.clear();
 gpat.seekg(0, gpat.beg);

 cerr << "INFO: read values to permute: " << data.size() << endl;

 srand (time(NULL));

 if(gpat.is_open()){

   while(getline(gpat, line)){
     vector<string> region = split(line, "\t");

     double value = atof(region[4].c_str());
     
     if(value < 0){
       value = 0;
     }
     

     double suc   = 0;
     double per   = 0;
     int    datas = data.size();
     double pv = (1.0 / globalOpts.npermutation);     

     while( suc < globalOpts.nsuc && per < globalOpts.npermutation){
       per += 1.0;
       
       int r = rand() % datas;

       if(value < data[r]){
	 suc += 1;
       }
     }
     if(suc > 0){
       pv = suc / per;
     }
     cout << line << "\t" << suc << "\t" << per << "\t" << pv << endl;
   }   
 }
 else{
   cerr << "FATAL: coult not open file: " << globalOpts.file << endl;
   exit(1);
 }


return 0;
}
