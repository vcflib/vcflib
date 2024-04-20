% VCFREMOVEABERRANTGENOTYPES(1) vcfremoveaberrantgenotypes (vcflib) | vcfremoveaberrantgenotypes (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfremoveaberrantgenotypes**

# SYNOPSIS

**vcfremoveaberrantgenotypes** <vcf file>

# DESCRIPTION

strips samples which are homozygous but have observations implying heterozygosity. Remove samples for which the reported genotype (GT) and observation counts disagree (AO, RO).



# OPTIONS

```


Type: transformation

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

[vcfremoveaberrantgenotypes.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfremoveaberrantgenotypes.cpp)

# LICENSE

Copyright 2011-2024 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2024 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
