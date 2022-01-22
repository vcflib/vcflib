% NORMALIZE-IHS(1) normalize-iHS (vcflib) | normalize-iHS (VCF genotype)
% Erik Garrison and vcflib contributors

# NAME

**normalize-iHS**

# SYNOPSIS

normalizeHS -s 0.01 -f input.txt

# DESCRIPTION

normalizes iHS or XP-EHH scores.



# OPTIONS

```




A cross-population extended haplotype homozygosity (XP-EHH) score is
directional: a positive score suggests selection is likely to have
happened in population A, whereas a negative score suggests the same
about population B. See for example
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2687721/


Output : **normalize-iHS** adds one additional column to input (normalized score).
required: -f            -- Output from iHS or XPEHH 
optional: -s            -- Max AF diff for window [0.01]

Type: genotype



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

[normalize-iHS.cpp](https://github.com/vcflib/vcflib/blob/master/src/normalize-iHS.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
