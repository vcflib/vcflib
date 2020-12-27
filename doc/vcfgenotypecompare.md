% VCFGENOTYPECOMPARE(1) vcfgenotypecompare (vcflib) | vcfgenotypecompare (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**vcfgenotypecompare**

# SYNOPSIS

**vcfgenotypecompare** <other-genotype-tag> <vcf file>

# DESCRIPTION

adds statistics to the INFO field of the vcf file describing the amount of discrepancy between the genotypes (GT) in the vcf file and the genotypes reported in the <other-genotype-tag>. use this after vcfannotategenotypes to get correspondence statistics for two vcfs.



# OPTIONS

```


Type: statistics

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcfgenotypecompare.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfgenotypecompare.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
