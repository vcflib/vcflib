#include "Variant.h"
#include "split.h"
#include "var.hpp"

#include <string>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

using namespace vcf;

int main(int argc, char** argv) {

  string filename = argv[1];

  VariantCallFile variantFile;

  variantFile.open(filename);

  vector<string> headerLines = split (variantFile.header, "\n");
  
  for(vector<string>::iterator it = headerLines.begin(); it != headerLines.end(); it++){

    //    cerr << "h:" <<  (*it) << endl;
    
  if((*it).substr(0,8) == "##contig"){
    string contigInfo = (*it).substr(13, (*it).length() -14);

    vector<string> info = split(contigInfo, ",");
    cout << info[0] << "\t" << info[1].substr(7) << endl;
    
  }


  }

}
