/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#ifndef gpatInfo_H
#define gpatInfo_H

#include <iostream>

void printBasicVersion(void){
  std::cout << "vcflib" << "\t" << VCFLIB_VERSION << std::endl;
}

void printVersion(void){

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "              This is a vcflib::GPAT++ tool           " << std::endl;

  std::cerr << "-Version       : " << VCFLIB_VERSION << std::endl;
  std::cerr << "-Contact       : zev.kronenberg [at] gmail.com " << std::endl;
  std::cerr << "-Notes         : If you find a bug, please open a report on github!" << std::endl;
  std::cerr << "-Support       : Please post questions to biostars.org             " << std::endl;
  std::cerr << "-Contribution  :                                                   " << std::endl;
  std::cerr << "                 Zev Kronenberg (UW Genome Sciences)               " << std::endl;
  std::cerr << "                 Mark Yandell   (UU Human genetics)                " << std::endl;
  std::cerr << "                 Mike Shapiro   (UU Biology)                       " << std::endl;
  std::cerr << "                 EJ Osborne     (UU Human genetics)                " << std::endl;
  std::cerr << "                 Brett Kennedy  (UU Human genetics)                " << std::endl;
  std::cerr << "                 Daniel Ence    (UU Human genetics)                " << std::endl;
  std::cerr << "                 Erik Garrison  (Wellcome Trust Sanger Institute)  " << std::endl;
  std::cerr << "                 Travis Collier (UC Davis)                         " << std::endl;
  std::cerr << "                 -     Your name goes here       -'                " << std::endl;

  std::cerr << "------------------------------------------------------" << std::endl;
}

#endif
