% VCFGLXGT(1) vcfglxgt (vcflib) | vcfglxgt (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfglxgt**

# SYNOPSIS

**vcfglxgt** [options] <vcf file>

# DESCRIPTION

Set genotypes using the maximum genotype likelihood for each sample.



# OPTIONS

```

options:
    -n, --fix-null-genotypes   only apply to null and partly-null genotypes



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

[vcfglxgt.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfglxgt.cpp)

# LICENSE

Copyright 2011-2025 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2025 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
