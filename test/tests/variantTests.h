#include "gtest/gtest.h"
#include "Variant.h"
#include <iostream>
#include <string>
#include <vector>

TEST(VariantCallFile, open){

  vcflib::VariantCallFile variantFile;

  std::string filename = "../samples/sample.vcf";

  variantFile.open(filename);
  
  ASSERT_TRUE(variantFile.is_open());


};


TEST(VariantCallFile, recordCountUncompressed){
  vcflib::VariantCallFile variantFile;
  
  std::string filename = "../samples/sample.vcf";
  variantFile.open(filename);
  vcflib::Variant var(variantFile);
  
  long int count = 0;

  while (variantFile.getNextVariant(var)) {
    count+= 1;
  }
  ASSERT_EQ(count, 9);
}

TEST(VariantCallFile, recordCountCompressed){
  vcflib::VariantCallFile variantFile;

  std::string filename = "../samples/sample.compressed.vcf.gz";
  variantFile.open(filename);
  vcflib::Variant var(variantFile);

  long int count = 0;
  while (variantFile.getNextVariant(var)) {
    count+= 1;
  }
  ASSERT_EQ(count, 9);
}

TEST(VariantCallFile, sampleNamesCompressed){
  vcflib::VariantCallFile variantFile;

  std::string filename = "../samples/sample.compressed.vcf.gz";
  variantFile.open(filename);

  int sampleSize = variantFile.sampleNames.size();

  ASSERT_EQ(sampleSize, 3);

  std::vector<std::string> names;
  names.push_back("NA00001");
  names.push_back("NA00002");
  names.push_back("NA00003");
  
  int counter = 0;

  for(std::vector<std::string>::iterator it 
	= variantFile.sampleNames.begin(); 
      it != variantFile.sampleNames.end(); it++){
    
    ASSERT_EQ(*it, names[counter]);
    
    counter+=1;

  }


  

}
