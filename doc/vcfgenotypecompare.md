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

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcfgenotypecompare.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfgenotypecompare.cpp)

# LICENSE

Copyright 2011-2025 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2025 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
