% VCFRANDOM(1) vcfrandom (vcflib) | vcfrandom (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

vcfrandom

# SYNOPSIS



# DESCRIPTION

##fileformat=VCFv4.0 ##source=vcfrandom ##reference=/d2/data/references/build_37/human_reference_v37.fa ##phasing=none ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data"> ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus"> ##INFO=<ID=AC,Number=1,Type=Integer,Description="Total number of alternate alleles in called genotypes"> ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes"> ##INFO=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the range (0,1]"> ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"> ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype"> ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth"> #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT bill one 1 . T G,T 100 . DP=60 GT:DP 0/1:52 one 2 . G T,A 100 . DP=14 GT:DP 0/1:86 one 3 . A G,C 100 . DP=87 GT:DP 0/1:45 one 4 . A G,A 100 . DP=79 GT:DP 0/1:98 one 5 . G A,A 100 . DP=28 GT:DP 0/1:51 one 6 . A G,C 100 . DP=68 GT:DP 0/1:98 one 7 . G G,A 100 . DP=24 GT:DP 0/1:29 one 8 . G A,C 100 . DP=78 GT:DP 0/1:86 one 9 . G T,T 100 . DP=28 GT:DP 0/1:19





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
