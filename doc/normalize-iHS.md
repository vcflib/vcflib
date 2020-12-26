% NORMALIZE-IHS(1) normalize-iHS (vcflib) | normalize-iHS (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

normalize-iHS

# SYNOPSIS

usage: normalizeHS -s 0.01 -f input.txt

# DESCRIPTION

normalizes iHS or XP-EHH scores Output : normalize-iHS adds one additional column to input (normalized score). required: -f -- Output from iHS or XPEHH optional: -s -- Max AF diff for window [0.01]

# OPTIONS

```



```



# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[normalize-iHS.cpp](https://github.com/vcflib/vcflib/blob/master/src/normalize-iHS.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
