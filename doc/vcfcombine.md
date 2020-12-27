% VCFCOMBINE(1) vcfcombine (vcflib) | vcfcombine (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

**vcfcombine**

# SYNOPSIS

**vcfcombine** [vcf file] [vcf file] ...

# DESCRIPTION

Combines VCF files positionally, combining samples when sites and alleles are identical. Any number of VCF files may be combined. The INFO field and other columns are taken from one of the files which are combined when records in multiple files match. Alleles must have identical ordering to be combined into one record. If they do not, multiple records will be emitted.



# OPTIONS

```


options:
    -h --help           This text.
    -v --version        Print version.
    -r --region REGION  A region specifier of the form chrN:x-y to bound the merge

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcfcombine.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfcombine.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
