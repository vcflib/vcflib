% VCFGLXGT(1) vcfglxgt (vcflib) | vcfglxgt (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

vcfglxgt

# SYNOPSIS

usage: ./build/vcfglxgt [options] <vcf file>

# DESCRIPTION

options: -n, --fix-null-genotypes only apply to null and partly-null genotypes

# OPTIONS

```


Set genotypes using the maximum genotype likelihood for each sample.

```

# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcfglxgt.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfglxgt.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
