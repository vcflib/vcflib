#include "gtest/gtest.h"
#include "Variant.h"
#include <iostream>

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
