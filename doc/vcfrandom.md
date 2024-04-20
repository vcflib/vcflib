% VCFRANDOM(1) vcfrandom (vcflib) | vcfrandom (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**vcfrandom**

# SYNOPSIS

**vcfrandom**

# DESCRIPTION

Generate a random VCF file





# EXAMPLES

```

Example:

    **vcfrandom**

##fileformat=VCFv4.0
##source=**vcfrandom**
##reference=/d2/data/references/build_37/human_reference_v37.fa
##phasing=none
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
##INFO=<ID=AC,Number=1,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the range (0,1]">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  bill
one     1       .       G       G,A     100     .       DP=83   GT:DP   0/1:1
one     2       .       G       G,A     100     .       DP=3    GT:DP   0/1:49
one     3       .       G       C,T     100     .       DP=5    GT:DP   0/1:12
one     4       .       C       G,T     100     .       DP=51   GT:DP   0/1:60
one     5       .       A       T,A     100     .       DP=31   GT:DP   0/1:89
one     6       .       T       T,A     100     .       DP=56   GT:DP   0/1:60
one     7       .       T       A,C     100     .       DP=78   GT:DP   0/1:75
one     8       .       T       G,A     100     .       DP=73   GT:DP   0/1:78
one     9       .       C       C,G     100     .       DP=42   GT:DP   0/1:67


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

[vcfrandom.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfrandom.cpp)

# LICENSE

Copyright 2011-2024 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2024 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
