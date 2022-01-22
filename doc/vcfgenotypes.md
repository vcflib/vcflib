% VCFGENOTYPES(1) vcfgenotypes (vcflib) | vcfgenotypes (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**vcfgenotypes**

# SYNOPSIS

**vcfgenotypes** <vcf file>

# DESCRIPTION

Report the genotypes for each sample, for each variant in the VCF. Convert the numerical represenation of genotypes provided by the GT field to a human-readable genotype format.



# OPTIONS

```




```





# EXAMPLES

```

Example:

      **vcfgenotypes** samples/sample.vcf

19      111     A       C       A,C     NA00001:A/A     NA00002:A/A     NA00003:A/C
19      112     A       G       A,G     NA00001:A/A     NA00002:A/A     NA00003:A/G
20      14370   G       A       G,A     NA00001:G/G     NA00002:G/A     NA00003:A/A
20      17330   T       A       T,A     NA00001:T/T     NA00002:T/A     NA00003:T/T
20      1110696 A       G,T     A,G,T   NA00001:G/T     NA00002:G/T     NA00003:T/T
20      1230237 T       .       T,.     NA00001:T/T     NA00002:T/T     NA00003:T/T
20      1234567 G       GA,GAC  G,GA,GAC        NA00001:G/GA    NA00002:G/GAC   NA00003:GA/GA
20      1235237 T       .       T,.     NA00001:T/T     NA00002:T/T     NA00003:./.
X       10      AC      A,ATG   AC,A,ATG        NA00001:AC      NA00002:AC/A    NA00003:AC/ATG

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

[vcfgenotypes.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfgenotypes.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
