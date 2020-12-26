% VCFRANDOM(1) vcfrandom (vcflib) | vcfrandom (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

vcfrandom

# SYNOPSIS



# DESCRIPTION

##fileformat=VCFv4.0 ##source=vcfrandom ##reference=/d2/data/references/build_37/human_reference_v37.fa ##phasing=none ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data"> ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus"> ##INFO=<ID=AC,Number=1,Type=Integer,Description="Total number of alternate alleles in called genotypes"> ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes"> ##INFO=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the range (0,1]"> ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"> ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype"> ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth"> #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT bill one 1 . A A,C 100 . DP=39 GT:DP 0/1:55 one 2 . C A,A 100 . DP=49 GT:DP 0/1:48 one 3 . G C,G 100 . DP=19 GT:DP 0/1:48 one 4 . T C,T 100 . DP=82 GT:DP 0/1:34 one 5 . A T,C 100 . DP=34 GT:DP 0/1:75 one 6 . C A,G 100 . DP=73 GT:DP 0/1:43 one 7 . G G,G 100 . DP=86 GT:DP 0/1:25 one 8 . G C,G 100 . DP=5 GT:DP 0/1:54 one 9 . C T,C 100 . DP=20 GT:DP 0/1:62

# OPTIONS

```



```

# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcfrandom.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfrandom.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
