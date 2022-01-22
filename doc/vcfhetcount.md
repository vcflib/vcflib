% VCFHETCOUNT(1) vcfhetcount (vcflib) | vcfhetcount (VCF metrics)
% Erik Garrison and vcflib contributors

# NAME

**vcfhetcount**

# SYNOPSIS

**vcfhetcount** <vcf file>

# DESCRIPTION

Calculate the heterozygosity rate: count the number of alternate alleles in heterozygous genotypes in all records in the vcf file



# OPTIONS

```

outputs a count for each individual in the file

Type: metrics

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

[vcfhetcount.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfhetcount.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
