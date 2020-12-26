% VCFCREATEMULTI(1) vcfcreatemulti (vcflib) | vcfcreatemulti (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

vcfcreatemulti

# SYNOPSIS

usage: ./build/vcfcreatemulti [options] [file]

# DESCRIPTION

If overlapping alleles are represented across multiple records, merge them into a single record. Currently only for indels.





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcfcreatemulti.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfcreatemulti.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
