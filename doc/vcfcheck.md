% VCFCHECK(1) vcfcheck (vcflib) | vcfcheck (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

**vcfcheck**

# SYNOPSIS

./build/**vcfcheck** [options] <vcf file>

# DESCRIPTION

options: -f, --fasta-reference FASTA reference file to use to obtain primer sequences -x, --exclude-failures If a record fails, don't print it. Otherwise do. -k, --keep-failures Print if the record fails, otherwise not. -h, --help Print this message. -v, --version Print version.



# OPTIONS

```


Verifies that the VCF REF field matches the reference as described.

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcfcheck.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfcheck.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
