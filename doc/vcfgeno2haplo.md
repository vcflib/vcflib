% VCFGENO2HAPLO(1) vcfgeno2haplo (vcflib) | vcfgeno2haplo (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfgeno2haplo**

# SYNOPSIS

**vcfgeno2haplo** [options] [<vcf file>]

# DESCRIPTION

Convert genotype-based phased alleles within --window-size into haplotype alleles. Will break haplotype construction when encountering non-phased genotypes on input.



# OPTIONS

```

options:
    -h, --help              Print this message
    -v, --version           Print version
    -r, --reference FILE    FASTA reference file
    -w, --window-size N     Merge variants at most this many bp apart (default 30)
    -o, --only-variants     Don't output the entire haplotype, just concatenate
                            REF/ALT strings (delimited by ":")



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

[vcfgeno2haplo.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfgeno2haplo.cpp)

# LICENSE

Copyright 2011-2021 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
