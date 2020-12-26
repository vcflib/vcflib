% VCFGLBOUND(1) vcfglbound (vcflib) | vcfglbound (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

vcfglbound

# SYNOPSIS

usage: ./build/vcfglbound [options] <vcf file>

# DESCRIPTION

options: -b, --bound N Bound GLs to this limit. -x, --exclude-broken If GLs are > 0, remove site.

# OPTIONS

```


Adjust GLs so that the maximum GL is 0 by dividing all GLs for each sample by the max.
Then cap (bound) at N (e.g. -10).

```

# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcfglbound.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfglbound.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
