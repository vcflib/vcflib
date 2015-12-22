#include "Variant.h"
#include "split.h"
#include "var.hpp"

#include <string>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

using namespace vcflib;

int main(int argc, char** argv) {

  string filename = argv[1];

  VariantCallFile variantFile;

  variantFile.open(filename);

  vector<string> headerLines = split (variantFile.header, "\n");
  
  for(vector<string>::iterator it = headerLines.begin(); it != headerLines.end(); it++){

    //    cerr << "h:" <<  (*it) << endl;
    
  if((*it).substr(0,8) == "##contig"){
    string contigInfo = (*it).substr(10, (*it).length() -11);
    //    cerr << contigInfo << endl;
    vector<string> info = split(contigInfo, ",");
    for(vector<string>::iterator sub = info.begin(); sub != info.end(); sub++){
      //      cerr << "s:" << (*sub) << endl;
      vector<string> subfield = split((*sub), "=");
      if(subfield[0] == "ID"){
	cout << subfield[1] << "\t";
      }
      if(subfield[0] == "length"){
	cout << subfield[1] << endl;
      }
    }
    
  }


  }

}
