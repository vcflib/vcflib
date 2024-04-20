% VCFGLBOUND(1) vcfglbound (vcflib) | vcfglbound (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfglbound**

# SYNOPSIS

**vcfglbound** [options] <vcf file>

# DESCRIPTION

Adjust GLs so that the maximum GL is 0 by dividing all GLs for each sample by the max.



# OPTIONS

```


Then cap (bound) at N (e.g. -10).options:
    -b, --bound N          Bound GLs to this limit.
    -x, --exclude-broken   If GLs are > 0, remove site.


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

[vcfglbound.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfglbound.cpp)

# LICENSE

Copyright 2011-2024 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2024 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
