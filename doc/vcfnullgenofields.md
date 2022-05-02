% VCFNULLGENOFIELDS(1) vcfnullgenofields (vcflib) | vcfnullgenofields (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfnullgenofields**

# SYNOPSIS

**vcfnullgenofields** [options] <vcf file>

# DESCRIPTION

Makes the FORMAT for each variant line the same (uses all the FORMAT fields described in the header). Fills out per-sample fields to match FORMAT. Expands GT values of '.' with number of alleles based on ploidy (eg: './.' for dipolid).



# OPTIONS

```


options:

    -p, --ploidy N   the polidy of missing/null GT fields (default=2)
    -L, --expand_GL  fill in missing GL fields with 0 values (eg: 0,0,0 for diploid 2 alleles)

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

[vcfnullgenofields.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfnullgenofields.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
