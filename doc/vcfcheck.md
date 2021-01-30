% VCFCHECK(1) vcfcheck (vcflib) | vcfcheck (VCF metrics)
% Erik Garrison and vcflib contributors

# NAME

**vcfcheck**

# SYNOPSIS

**vcfcheck** [options] <vcf file>

# DESCRIPTION

Validate integrity and identity of the VCF by verifying that the VCF record's REF matches a given reference file.



# OPTIONS

```

options:
    -f, --fasta-reference  FASTA reference file to use to obtain primer sequences
    -x, --exclude-failures If a record fails, don't print it.  Otherwise do.
    -k, --keep-failures    Print if the record fails, otherwise not.
    -h, --help       Print this message.
    -v, --version    Print version.


Type: metrics

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

[vcfcheck.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfcheck.cpp)

# LICENSE

Copyright 2011-2021 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
