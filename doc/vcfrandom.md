% VCFRANDOM(1) vcfrandom (vcflib) | vcfrandom (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

vcfrandom

# SYNOPSIS



# DESCRIPTION

##fileformat=VCFv4.0 ##source=vcfrandom ##reference=/d2/data/references/build_37/human_reference_v37.fa ##phasing=none ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data"> ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus"> ##INFO=<ID=AC,Number=1,Type=Integer,Description="Total number of alternate alleles in called genotypes"> ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes"> ##INFO=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the range (0,1]"> ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"> ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype"> ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth"> #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT bill one 1 . C T,G 100 . DP=24 GT:DP 0/1:28 one 2 . T A,A 100 . DP=78 GT:DP 0/1:15 one 3 . G G,G 100 . DP=63 GT:DP 0/1:25 one 4 . A G,G 100 . DP=40 GT:DP 0/1:84 one 5 . C C,T 100 . DP=23 GT:DP 0/1:14 one 6 . T C,C 100 . DP=17 GT:DP 0/1:12 one 7 . G A,C 100 . DP=78 GT:DP 0/1:28 one 8 . T C,T 100 . DP=80 GT:DP 0/1:61 one 9 . G A,A 100 . DP=15 GT:DP 0/1:5

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
