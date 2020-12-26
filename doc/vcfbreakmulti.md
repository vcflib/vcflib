% VCFBREAKMULTI(1) vcfbreakmulti (vcflib) | vcfbreakmulti (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

vcfbreakmulti

# SYNOPSIS

usage: ./build/vcfbreakmulti [options] [file]

# DESCRIPTION

If multiple alleles are specified in a single record, break the record into multiple lines, preserving allele-specific INFO fields.

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

[vcfbreakmulti.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfbreakmulti.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
